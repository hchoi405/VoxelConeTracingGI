#include "WrapBorderPass.h"
#include "engine/rendering/architecture/RenderPipeline.h"
#include "engine/resource/ResourceManager.h"
#include "engine/rendering/Texture3D.h"
#include "Globals.h"

WrapBorderPass::WrapBorderPass()
    : RenderPass("WrapBorderPass")
{
    m_shader = ResourceManager::getComputeShader("shaders/voxelConeTracing/copyWrappedBorder.comp");
}

void WrapBorderPass::update()
{
    // Fetch the data
    auto voxelOpacity = m_renderPipeline->fetchPtr<Texture3D>("VoxelOpacity");
    auto voxelRadiance = m_renderPipeline->fetchPtr<Texture3D>("VoxelRadiance");

    copyBorder(voxelOpacity, false);
    copyBorder(voxelRadiance, false);

    auto virtualVoxelOpacity = m_renderPipeline->fetchPtr<Texture3D>("VirtualVoxelOpacity");
    auto virtualVoxelRadiance = m_renderPipeline->fetchPtr<Texture3D>("VirtualVoxelRadiance");

    copyBorder(virtualVoxelOpacity, true);
    copyBorder(virtualVoxelRadiance, true);
}

void WrapBorderPass::copyBorder(Texture3D* texture, bool isVirtual) const
{
    m_shader->bind();

    int voxelResolution = (isVirtual)? VIRTUAL_VOXEL_RESOLUTION : VOXEL_RESOLUTION;
    int clipRegionCount = (isVirtual)? VIRTUAL_CLIP_REGION_COUNT : CLIP_REGION_COUNT;

    m_shader->setInt("u_clipmapResolution", voxelResolution);
    m_shader->setInt("u_clipmapResolutionWithBorder", voxelResolution + 2);
    m_shader->setInt("u_faceCount", FACE_COUNT);
    m_shader->setInt("u_clipmapCount", clipRegionCount);
    m_shader->bindImage3D(*texture, "u_image", GL_READ_WRITE, GL_RGBA8, 0);

    float borderWidth2 = 2.0f;
    GLuint groupCount = GLuint(ceil((voxelResolution + borderWidth2) / 8.0f));
    m_shader->dispatchCompute(groupCount, groupCount, groupCount);
    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
}
