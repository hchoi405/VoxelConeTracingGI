#include "MaterialEstimationPass.h"
#include <GL/glew.h>
#include <engine/camera/CameraComponent.h>
#include <engine/rendering/Texture3D.h>
#include <engine/rendering/architecture/RenderPipeline.h>
#include <engine/resource/ResourceManager.h>
#include <engine/ecs/ECS.h>
#include "../voxelConeTracing/settings/VoxelConeTracingSettings.h"
#include "../voxelConeTracing/voxelization.h"
#include "../voxelConeTracing/Downsampler.h"
#include "engine/util/QueryManager.h"
#include "../voxelConeTracing/VoxelConeTracing.h"
#include "engine/rendering/util/ImageCleaner.h"
#include "engine/util/ECSUtil/ECSUtil.h"
#include "../voxelConeTracing/ClipmapUpdatePolicy.h"

MaterialEstimationPass::MaterialEstimationPass()
    : RenderPass("MaterialEstimationPass")
{
    // m_materialEstimaionShader = ResourceManager::getShader("shaders/tvcg17/materialEstimation.vert",
    //                                                        "shaders/tvcg17/materialEstimation.frag", "shaders/tvcg17/materialEstimation.geom");
    m_materialEstimaionShader = ResourceManager::getComputeShader("shaders/tvcg17/materialEstimation.comp");
}

void MaterialEstimationPass::update()
{
    auto voxelRadiance = m_renderPipeline->fetchPtr<Texture3D>("VoxelRadiance");
    auto voxelOpacity = m_renderPipeline->fetchPtr<Texture3D>("VoxelOpacity");
    auto voxelNormal = m_renderPipeline->fetchPtr<Texture3D>("VoxelNormal");
    auto voxelReflectance = m_renderPipeline->fetchPtr<Texture3D>("VoxelReflectance");
    auto clipRegions = m_renderPipeline->fetchPtr<std::vector<VoxelRegion>>("ClipRegions");

    if (m_initializing) {
        estimateMaterial(voxelRadiance, voxelOpacity, voxelNormal, voxelReflectance, clipRegions->at(0), VOXEL_RESOLUTION);
        // downsample(voxelRadiance, *clipRegions, CLIP_REGION_COUNT, VOXEL_RESOLUTION);
        // downsample(voxelNormal, *clipRegions, CLIP_REGION_COUNT, VOXEL_RESOLUTION);
    }

#ifdef VIRTUAL
    
#endif

    // m_initializing = false;
}

void MaterialEstimationPass::estimateMaterial(Texture3D *voxelRadiance, Texture3D *voxelOpacity, Texture3D *voxelNormal, Texture3D *voxelReflectance,
                                              VoxelRegion regionL0, int voxelResolution) const
{
    m_materialEstimaionShader->bind();
    m_materialEstimaionShader->bindTexture3D(*voxelRadiance, "u_voxelRadiance", 0);
    m_materialEstimaionShader->bindTexture3D(*voxelOpacity, "u_voxelOpacity", 1);
    m_materialEstimaionShader->bindTexture3D(*voxelNormal, "u_voxelNormal", 2);
    m_materialEstimaionShader->bindImage3D(*voxelReflectance, "u_voxelReflectance", GL_READ_WRITE, GL_RGBA8, 3);

    auto camera = m_renderPipeline->getCamera();
    m_materialEstimaionShader->setMatrix("u_viewProjInv", camera->viewProjInv());
    m_materialEstimaionShader->setMatrix("u_viewProj", camera->viewProj());
    auto dasan = ECS::getEntityByName("dasan613.obj");
    m_materialEstimaionShader->setMatrix("u_model", dasan.getComponent<Transform>()->getLocalToWorldMatrix());

    m_materialEstimaionShader->setFloat("u_traceStartOffset", GI_SETTINGS.traceStartOffset);
    m_materialEstimaionShader->setFloat("u_stepFactor", GI_SETTINGS.stepFactor);
    m_materialEstimaionShader->setVectori("u_imageMin", regionL0.getMinPosImage(regionL0.extent));
    m_materialEstimaionShader->setVectori("u_regionMin", regionL0.minPos);
    // std::cout << "u_imageMin: " << regionL0.getMinPosImage(regionL0.extent) << std::endl;
    // std::cout << "u_regionMin: " << regionL0.minPos << std::endl;
    // std::cout << "getMaxPos(): " << regionL0.getMaxPos() << std::endl;
    m_materialEstimaionShader->setFloat("u_voxelSizeL0", regionL0.voxelSize);
    // std::cout << "u_voxelSizeL0: " << regionL0.voxelSize << std::endl;
    // std::cout << "u_worldMin from voxel: " << glm::vec3(regionL0.minPos) * regionL0.voxelSize << std::endl;
    // std::cout << "u_worldMin from voxel: " << glm::vec3(regionL0.getMaxPos()) * regionL0.voxelSize << std::endl;
    m_materialEstimaionShader->setVector("u_volumeCenterL0", regionL0.getCenterPosWorld());
    m_materialEstimaionShader->setFloat("u_occlusionDecay", GI_SETTINGS.occlusionDecay);
    m_materialEstimaionShader->setFloat("u_indirectDiffuseIntensity", GI_SETTINGS.indirectDiffuseIntensity);
    m_materialEstimaionShader->setUnsignedInt("u_irradianceOnly", DEBUG_SETTINGS.irradianceOnly.value? 1 : 0);

    m_materialEstimaionShader->setInt("u_clipmapResolution", voxelResolution);
    m_materialEstimaionShader->setInt("u_clipmapResolutionWithBorder", voxelResolution + 2);
    m_materialEstimaionShader->setInt("u_clipmapLevel", 0);

    QueryManager::beginElapsedTime(QueryTarget::GPU, "Estimate Material");

    GLuint groupCount = voxelResolution / 8;
    m_materialEstimaionShader->dispatchCompute(groupCount, groupCount, groupCount);

    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

    QueryManager::endElapsedTime(QueryTarget::GPU, "Estimate Material");
}

void MaterialEstimationPass::downsample(Texture3D *voxelRadiance, std::vector<VoxelRegion> &clipRegions, int clipRegionCount, int voxelResolution) const
{
    QueryManager::beginElapsedTime(QueryTarget::GPU, "Reflectance Downsampling");

    Downsampler::downsample(voxelRadiance, &clipRegions, clipRegionCount, voxelResolution);
    
    QueryManager::endElapsedTime(QueryTarget::GPU, "Reflectance Downsampling");
}