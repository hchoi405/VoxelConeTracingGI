#pragma once
#include <memory>
#include <glm/glm.hpp>
#include <vector>
#include "VoxelRegion.h"
#include "ClipmapUpdatePolicy.h"

class Texture3D;
class Shader;

class CopyAlpha
{
public:
    static void init();

    static void copyAlpha(Texture3D *voxelRadiance, Texture3D *voxelOpacity, int clipLevel, int voxelResolution);
    static void copyAlpha(Texture3D *voxelRadiance, Texture3D *voxelOpacity, int voxelResolution, ClipmapUpdatePolicy *clipmapUpdatePolicy);
    

private:
    static std::shared_ptr<Shader> m_copyAlphaShader;
};
