#pragma once
#include <engine/rendering/shader/Shader.h>
#include <engine/rendering/geometry/Mesh.h>
#include <engine/rendering/architecture/RenderPass.h>
#include "../voxelConeTracing/voxelization.h"
#include "../voxelConeTracing/Globals.h"
#include "../voxelConeTracing/ClipmapUpdatePolicy.h"

class Framebuffer;
class Texture3D;
class BBox;

class MaterialEstimationPass : public RenderPass
{
public:
    explicit MaterialEstimationPass();

    void update() override;

private:
    void estimateMaterial(Texture3D *voxelRadiance, Texture3D *voxelOpacity, Texture3D *voxelNormal, Texture3D *voxelReflectance,
                          VoxelRegion regionL0, int voxelResolution) const;
    void estimateVirtualMaterial(Texture3D *virtualVoxelOpacity, Texture3D *voxelDiffuse, Texture3D *voxelSpecularA,
                                 VoxelRegion regionL0, int virtualVoxelResolution) const;

    void downsample(Texture3D *voxelRadiance, std::vector<VoxelRegion> &clipRegions, int clipRegionCount, int voxelResolution) const;

private:
    std::shared_ptr<Shader> m_materialEstimaionShader;
    std::shared_ptr<Shader> m_virutalMaterialVoxelizeShader;
    VoxelizationMode m_voxelizationMode{VoxelizationMode::CONSERVATIVE};

    bool m_initializing{true};
};
