#include "RadianceInjectionPass.h"
#include <GL/glew.h>
#include <engine/camera/CameraComponent.h>
#include <engine/rendering/Texture3D.h>
#include <engine/rendering/architecture/RenderPipeline.h>
#include <engine/resource/ResourceManager.h>
#include <engine/ecs/ECS.h>
#include "settings/VoxelConeTracingSettings.h"
#include "voxelization.h"
#include "Downsampler.h"
#include "engine/util/QueryManager.h"
#include "VoxelConeTracing.h"
#include "engine/rendering/util/ImageCleaner.h"
#include "engine/util/ECSUtil/ECSUtil.h"
#include "ClipmapUpdatePolicy.h"
#include "CopyAlpha.h"

RadianceInjectionPass::RadianceInjectionPass()
    : RenderPass("RadianceInjectionPass")
{
    m_conservativeVoxelizationShader = ResourceManager::getShader("shaders/voxelConeTracing/injectLightByConservativeVoxelization.vert",
                                                                  "shaders/voxelConeTracing/injectLightByConservativeVoxelization.frag", "shaders/voxelConeTracing/injectLightByConservativeVoxelization.geom");

    m_msaaVoxelizationShader = ResourceManager::getShader("shaders/voxelConeTracing/injectLightByMSAAVoxelization.vert",
                                                          "shaders/voxelConeTracing/injectLightByMSAAVoxelization.frag", "shaders/voxelConeTracing/injectLightByMSAAVoxelization.geom");

    m_pointCloudVoxelizationShader = ResourceManager::getShader("shaders/voxelConeTracing/injectLightByPointCloud.vert",
                                                                "shaders/voxelConeTracing/injectLightByPointCloud.frag", "shaders/voxelConeTracing/injectLightByPointCloud.geom");

    m_cachedClipRegions.resize(CLIP_REGION_COUNT);
    m_virtualCachedClipRegions.resize(VIRTUAL_CLIP_REGION_COUNT);
}

void RadianceInjectionPass::update()
{
    auto voxelRadiance = m_renderPipeline->fetchPtr<Texture3D>("VoxelRadiance");
    auto voxelOpacity = m_renderPipeline->fetchPtr<Texture3D>("VoxelOpacity");
    auto voxelNormal = m_renderPipeline->fetchPtr<Texture3D>("VoxelNormal");
    auto voxelReflectance = m_renderPipeline->fetchPtr<Texture3D>("VoxelReflectance");
    m_clipmapUpdatePolicy = m_renderPipeline->fetchPtr<ClipmapUpdatePolicy>("ClipmapUpdatePolicy");

    auto clipRegions = m_renderPipeline->fetchPtr<std::vector<VoxelRegion>>("ClipRegions");

    if (m_initializing)
    {
        m_cachedClipRegions = *clipRegions;
    }
    else
    {
        auto &levelsToUpdate = m_clipmapUpdatePolicy->getLevelsScheduledForUpdate();
        for (auto level : levelsToUpdate)
        {
            m_cachedClipRegions[level] = clipRegions->at(level);
        }
    }

    auto shader = getSelectedShader();
    shader->setInt("u_normalOnly", 0);
    injectByVoxelization(shader, voxelRadiance, voxelNormal, voxelReflectance, nullptr, m_voxelizationMode, *m_clipmapUpdatePolicy, VOXEL_RESOLUTION, m_cachedClipRegions);

    CopyAlpha::copyAlpha(voxelRadiance, voxelOpacity, VOXEL_RESOLUTION, m_clipmapUpdatePolicy);
    downsample(voxelRadiance, m_cachedClipRegions, CLIP_REGION_COUNT, VOXEL_RESOLUTION, m_clipmapUpdatePolicy);
    CopyAlpha::copyAlpha(voxelNormal, voxelOpacity, VOXEL_RESOLUTION, m_clipmapUpdatePolicy);
    downsample(voxelNormal, m_cachedClipRegions, CLIP_REGION_COUNT, VOXEL_RESOLUTION, m_clipmapUpdatePolicy);
    CopyAlpha::copyAlpha(voxelReflectance, voxelOpacity, VOXEL_RESOLUTION, m_clipmapUpdatePolicy);
    downsample(voxelReflectance, m_cachedClipRegions, CLIP_REGION_COUNT, VOXEL_RESOLUTION, m_clipmapUpdatePolicy);

#ifdef VIRTUAL
    auto virtualVoxelRadiance = m_renderPipeline->fetchPtr<Texture3D>("VirtualVoxelRadiance");
    auto virtualVoxelOpacity = m_renderPipeline->fetchPtr<Texture3D>("VirtualVoxelOpacity");
    auto virtualVoxelNormal = m_renderPipeline->fetchPtr<Texture3D>("VirtualVoxelNormal");
    auto virtualVoxelDiffuse = m_renderPipeline->fetchPtr<Texture3D>("VirtualVoxelDiffuse");
    auto virtualVoxelSpecularA = m_renderPipeline->fetchPtr<Texture3D>("VirtualVoxelSpecularA");
    m_virtualClipmapUpdatePolicy = m_renderPipeline->fetchPtr<ClipmapUpdatePolicy>("VirtualClipmapUpdatePolicy");

    auto virtualClipRegions = m_renderPipeline->fetchPtr<std::vector<VoxelRegion>>("VirtualClipRegions");

    if (m_initializing)
    {
        m_virtualCachedClipRegions = *virtualClipRegions;
    }
    else
    {
        auto &virtualLevelsToUpdate = m_virtualClipmapUpdatePolicy->getLevelsScheduledForUpdate();
        for (auto level : virtualLevelsToUpdate)
        {
            m_virtualCachedClipRegions[level] = virtualClipRegions->at(level);
        }
    }
    shader = m_conservativeVoxelizationShader.get();
    injectByVoxelization(shader, virtualVoxelRadiance, virtualVoxelNormal, virtualVoxelDiffuse, virtualVoxelSpecularA, VoxelizationMode::CONSERVATIVE, *m_virtualClipmapUpdatePolicy, VIRTUAL_VOXEL_RESOLUTION, m_virtualCachedClipRegions);
    CopyAlpha::copyAlpha(virtualVoxelRadiance, virtualVoxelOpacity, VIRTUAL_VOXEL_RESOLUTION, m_virtualClipmapUpdatePolicy);
    downsample(virtualVoxelRadiance, m_virtualCachedClipRegions, VIRTUAL_CLIP_REGION_COUNT, VIRTUAL_VOXEL_RESOLUTION, m_virtualClipmapUpdatePolicy);
    CopyAlpha::copyAlpha(virtualVoxelNormal, virtualVoxelOpacity, VIRTUAL_VOXEL_RESOLUTION, m_virtualClipmapUpdatePolicy);
    downsample(virtualVoxelNormal, m_virtualCachedClipRegions, VIRTUAL_CLIP_REGION_COUNT, VIRTUAL_VOXEL_RESOLUTION, m_virtualClipmapUpdatePolicy);
    CopyAlpha::copyAlpha(virtualVoxelDiffuse, virtualVoxelOpacity, VIRTUAL_VOXEL_RESOLUTION, m_virtualClipmapUpdatePolicy);
    downsample(virtualVoxelDiffuse, m_virtualCachedClipRegions, VIRTUAL_CLIP_REGION_COUNT, VIRTUAL_VOXEL_RESOLUTION, m_virtualClipmapUpdatePolicy);
    CopyAlpha::copyAlpha(virtualVoxelSpecularA, virtualVoxelOpacity, VIRTUAL_VOXEL_RESOLUTION, m_virtualClipmapUpdatePolicy);
    downsample(virtualVoxelSpecularA, m_virtualCachedClipRegions, VIRTUAL_CLIP_REGION_COUNT, VIRTUAL_VOXEL_RESOLUTION, m_virtualClipmapUpdatePolicy);
#endif

    m_initializing = false;
}

void RadianceInjectionPass::injectByVoxelization(Shader *shader, Texture3D *voxelRadiance, Texture3D *voxelNormal, Texture3D *voxelDiffuse,
                                                 Texture3D *voxelSpecularA, VoxelizationMode voxelizationMode, ClipmapUpdatePolicy &clipmapUpdatePolicy, int voxelResolution,
                                                 std::vector<VoxelRegion> &cachedClipRegions)
{
    static unsigned char zero[]{0, 0, 0, 0};

    QueryManager::beginElapsedTime(QueryTarget::GPU, "Clear Radiance Voxels");

    auto &levelsToUpdate = clipmapUpdatePolicy.getLevelsScheduledForUpdate();

    if (clipmapUpdatePolicy.getType() == ClipmapUpdatePolicy::Type::ALL_PER_FRAME)
    {
        glClearTexImage(*voxelRadiance, 0, GL_RGBA, GL_UNSIGNED_BYTE, zero);
    }
    else
    {
        for (auto level : levelsToUpdate)
        {
            auto clipLevel = static_cast<GLuint>(level);
            ImageCleaner::clear6FacesImage3D(*voxelRadiance, GL_RGBA8, glm::ivec3(0), glm::ivec3(voxelResolution), voxelResolution, clipLevel, 1);
        }

        glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
    }

    QueryManager::endElapsedTime(QueryTarget::GPU, "Clear Radiance Voxels");

    QueryManager::beginElapsedTime(QueryTarget::GPU, "Radiance Voxelization");
    VoxelizationDesc desc;
    desc.mode = voxelizationMode;
    desc.clipRegions = cachedClipRegions;
    desc.voxelizationShader = shader;
    if (voxelResolution == VIRTUAL_VOXEL_RESOLUTION) {
        auto virtualEntity = ECS::getEntityByName("virtualObject");
        desc.target = VoxelizationTarget::ENTITIES;
        desc.entities.push_back(virtualEntity);
    }
    else if (voxelResolution == VOXEL_RESOLUTION) {
        desc.target = VoxelizationTarget::POINTCLOUD;
    }
    desc.downsampleTransitionRegionSize = GI_SETTINGS.downsampleTransitionRegionSize;
    desc.voxelResolution = voxelResolution;
    Voxelizer *voxelizer = nullptr;
    if (voxelResolution == VIRTUAL_VOXEL_RESOLUTION)
        voxelizer = VoxelConeTracing::virtualVoxelizer();
    else if (voxelResolution == VOXEL_RESOLUTION)
        voxelizer = VoxelConeTracing::voxelizer();
    voxelizer->beginVoxelization(desc);
    shader->bindImage3D(*voxelRadiance, "u_voxelRadiance", GL_READ_WRITE, GL_R32UI, 0);
    shader->bindImage3D(*voxelNormal, "u_voxelNormal", GL_READ_WRITE, GL_R32UI, 1);
    shader->bindImage3D(*voxelDiffuse, "u_voxelDiffuse", GL_READ_WRITE, GL_R32UI, 2);
    if (voxelSpecularA) {
        shader->bindImage3D(*voxelSpecularA, "u_voxelSpecularA", GL_READ_WRITE, GL_R32UI, 3);
        shader->setUnsignedInt("u_noSpecular", 0);
    }
    else
        shader->setUnsignedInt("u_noSpecular", 1);

    for (auto level : levelsToUpdate)
    {
        voxelizer->voxelize(cachedClipRegions.at(level), level);
    }

    if (voxelResolution == VIRTUAL_VOXEL_RESOLUTION) {
        auto virtualEntity = ECS::getEntityByName("virtualObject");
        auto childTransforms = virtualEntity.getComponent<Transform>()->getChildren();
        for (auto childT: childTransforms) {
            childT->getOwner().getComponent<MeshRenderer>()->getMesh()->setRenderMode(GL_TRIANGLES, 0);
        }
    }

    voxelizer->endVoxelization(m_renderPipeline->getCamera()->getViewport());
    QueryManager::endElapsedTime(QueryTarget::GPU, "Radiance Voxelization");
}

void RadianceInjectionPass::downsample(Texture3D *voxelRadiance, std::vector<VoxelRegion> &cachedClipRegions, int clipRegionCount,
                                       int voxelResolution, ClipmapUpdatePolicy *clipmapUpdatePolicy) const
{
    QueryManager::beginElapsedTime(QueryTarget::GPU, "Radiance Downsampling");
    auto &levelsToUpdate = clipmapUpdatePolicy->getLevelsScheduledForUpdate();

    if (clipmapUpdatePolicy->getType() == ClipmapUpdatePolicy::Type::ALL_PER_FRAME)
    {
        Downsampler::downsample(voxelRadiance, &cachedClipRegions, clipRegionCount, voxelResolution);
    }
    else
    {
        for (auto level : levelsToUpdate)
        {
            if (level > 0)
                Downsampler::downsample(voxelRadiance, &cachedClipRegions, level, clipRegionCount, voxelResolution);
        }
    }

    QueryManager::endElapsedTime(QueryTarget::GPU, "Radiance Downsampling");
}

Shader *RadianceInjectionPass::getSelectedShader()
{
    switch (GI_SETTINGS.radianceInjectionMode)
    {
    case 0:
        m_voxelizationMode = VoxelizationMode::CONSERVATIVE;
        return m_conservativeVoxelizationShader.get();
    case 1:
        m_voxelizationMode = VoxelizationMode::MSAA;
        return m_msaaVoxelizationShader.get();
    case 2:
        // Use conservative voxelization to iterate all voxels
        // but inject radiance using point cloud
        m_voxelizationMode = VoxelizationMode::CONSERVATIVE;
        return m_pointCloudVoxelizationShader.get();
    default:
        assert(false);
        return nullptr;
    }
}
