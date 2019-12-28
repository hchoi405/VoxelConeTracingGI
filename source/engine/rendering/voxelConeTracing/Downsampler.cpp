#include "Downsampler.h"
#include "engine/resource/ResourceManager.h"
#include "engine/rendering/Texture3D.h"
#include "Globals.h"
#include "settings/VoxelConeTracingSettings.h"

std::shared_ptr<Shader> Downsampler::m_downsampleOpacityShader;
std::shared_ptr<Shader> Downsampler::m_downsampleShader;

void Downsampler::init()
{
    m_downsampleOpacityShader = ResourceManager::getComputeShader("shaders/voxelConeTracing/downsampleOpacity.comp");
    m_downsampleShader = ResourceManager::getComputeShader("shaders/voxelConeTracing/downsample3DImage.comp");
}

void Downsampler::downsampleOpacity(Texture3D* texture, const std::vector<VoxelRegion>* clipRegions, int clipmapLevel, int clipRegionCount, int voxelResolution)
{
    assert(clipmapLevel > 0 && clipmapLevel < clipRegionCount);

    m_downsampleOpacityShader->bind();

    m_downsampleOpacityShader->bindImage3D(*texture, "u_image", GL_READ_WRITE, GL_RGBA8, 0);

    m_downsampleOpacityShader->setInt("u_downsampleTransitionRegionSize", GI_SETTINGS.downsampleTransitionRegionSize);
    m_downsampleOpacityShader->setVectori("u_prevRegionMin", clipRegions->at(clipmapLevel - 1).minPos);
    m_downsampleOpacityShader->setInt("u_clipmapLevel", clipmapLevel);
    m_downsampleOpacityShader->setInt("u_clipmapResolution", voxelResolution);
    GLuint groupCount = GLuint(voxelResolution / 2 / 8);
    m_downsampleOpacityShader->dispatchCompute(groupCount, groupCount, groupCount);
    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
}

void Downsampler::downsample(Texture3D* image, const std::vector<VoxelRegion>* clipRegions, int clipRegionCount, int voxelResolution)
{
    m_downsampleShader->bind();
    m_downsampleShader->bindImage3D(*image, "u_image", GL_READ_WRITE, GL_RGBA8, 0);

    for (int i = 1; i < clipRegionCount; ++i)
    {
        m_downsampleShader->setInt("u_downsampleTransitionRegionSize", GI_SETTINGS.downsampleTransitionRegionSize);
        m_downsampleShader->setVectori("u_prevRegionMin", clipRegions->at(i - 1).minPos);
        m_downsampleShader->setInt("u_clipmapLevel", i);
        m_downsampleShader->setInt("u_clipmapResolution", voxelResolution);
        GLuint groupCount = GLuint(voxelResolution / 2 / 8);
        m_downsampleShader->dispatchCompute(groupCount, groupCount, groupCount);
        glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
    }
}

void Downsampler::downsample(Texture3D* image, const std::vector<VoxelRegion>* clipRegions, int clipmapLevel, int clipRegionCount, int voxelResolution)
{
    assert(clipmapLevel > 0 && clipmapLevel < clipRegionCount);

    m_downsampleShader->bind();
    m_downsampleShader->bindImage3D(*image, "u_image", GL_READ_WRITE, GL_RGBA8, 0);

    m_downsampleShader->setInt("u_downsampleTransitionRegionSize", GI_SETTINGS.downsampleTransitionRegionSize);
    m_downsampleShader->setVectori("u_prevRegionMin", clipRegions->at(clipmapLevel - 1).minPos);
    m_downsampleShader->setInt("u_clipmapLevel", clipmapLevel);
    m_downsampleShader->setInt("u_clipmapResolution", voxelResolution);
    GLuint groupCount = GLuint(voxelResolution / 2 / 8);
    m_downsampleShader->dispatchCompute(groupCount, groupCount, groupCount);
    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
}
