#include "CopyAlpha.h"
#include "engine/resource/ResourceManager.h"
#include "engine/rendering/Texture3D.h"
#include "Globals.h"
#include "settings/VoxelConeTracingSettings.h"

std::shared_ptr<Shader> CopyAlpha::m_copyAlphaShader;

void CopyAlpha::init()
{
    m_copyAlphaShader = ResourceManager::getComputeShader("shaders/voxelConeTracing/copyAlpha6Faces.comp");
}

void CopyAlpha::copyAlpha(Texture3D *voxelRadiance, Texture3D *voxelOpacity, int clipLevel, int voxelResolution)
{
    m_copyAlphaShader->bind();
    m_copyAlphaShader->bindTexture3D(*voxelOpacity, "u_srcTexture", 0);
    m_copyAlphaShader->bindImage3D(*voxelRadiance, "u_dstImage", GL_READ_WRITE, GL_RGBA8, 0);

    m_copyAlphaShader->setInt("u_clipmapResolution", voxelResolution);
    m_copyAlphaShader->setInt("u_clipmapResolutionWithBorder", voxelResolution + 2);
    m_copyAlphaShader->setInt("u_clipmapLevel", clipLevel);

    GLuint groupCount = voxelResolution / 8;
    m_copyAlphaShader->dispatchCompute(groupCount, groupCount, groupCount);
}

void CopyAlpha::copyAlpha(Texture3D *voxelRadiance, Texture3D *voxelOpacity, int voxelResolution, ClipmapUpdatePolicy *clipmapUpdatePolicy)
{
    auto &levelsToUpdate = clipmapUpdatePolicy->getLevelsScheduledForUpdate();

    // Copy alpha from opacity texture (this allows us to access just one texture during GI pass)
    for (auto level : levelsToUpdate)
    {
        copyAlpha(voxelRadiance, voxelOpacity, level, voxelResolution);
    }

    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
}