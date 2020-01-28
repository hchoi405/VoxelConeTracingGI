#pragma once
#include <engine/rendering/shader/Shader.h>
#include <engine/rendering/geometry/Mesh.h>
#include <engine/rendering/architecture/RenderPass.h>
#include "voxelization.h"
#include "Globals.h"
#include "ClipmapUpdatePolicy.h"

class Framebuffer;
class Texture3D;
class BBox;

struct DirLight
{
    glm::vec3 position;
    glm::vec3 direction;

    DirLight(const glm::vec3 &pos, const glm::vec3 &dir)
        : position(pos), direction(dir) {}

    DirLight() {}
};

class RadianceInjectionPass : public RenderPass
{
public:
    explicit RadianceInjectionPass();

    void update() override;
    Shader *getSelectedShader();

private:
    void injectByVoxelization(Shader *shader, Texture3D *voxelRadiance, Texture3D *voxelNormal, Texture3D *voxelDiffuse,
                              Texture3D *voxelSpecularA, VoxelizationMode voxelizationMode, ClipmapUpdatePolicy &clipmapUpdatePolicy, int voxelResolution,
                              std::vector<VoxelRegion> &cachedClipRegions);

    void downsample(Texture3D *voxelRadiance, std::vector<VoxelRegion> &cachedClipRegions, int clipRegionCount, int voxelResolution,
                    ClipmapUpdatePolicy *clipmapUpdatePolicy) const;

private:
    std::shared_ptr<Shader> m_conservativeVoxelizationShader;
    std::shared_ptr<Shader> m_msaaVoxelizationShader;
    std::shared_ptr<Shader> m_pointCloudVoxelizationShader;
    VoxelizationMode m_voxelizationMode{VoxelizationMode::CONSERVATIVE};

    ClipmapUpdatePolicy *m_clipmapUpdatePolicy{nullptr};
    ClipmapUpdatePolicy *m_virtualClipmapUpdatePolicy{nullptr};
    std::vector<VoxelRegion> m_cachedClipRegions;
    std::vector<VoxelRegion> m_virtualCachedClipRegions;
    bool m_initializing{true};
};
