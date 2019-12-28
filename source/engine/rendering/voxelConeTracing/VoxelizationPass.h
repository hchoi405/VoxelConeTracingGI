#pragma once
#include <engine/camera/CameraComponent.h>
#include <engine/rendering/shader/Shader.h>
#include <engine/rendering/architecture/RenderPass.h>
#include "Globals.h"
#include "engine/event/event.h"
#include "engine/event/EntityDeactivatedEvent.h"
#include "engine/event/EntityActivatedEvent.h"
#include "VoxelRegion.h"

class MeshRenderer;
class Texture3D;
class BBox;

class VoxelizationPass : public RenderPass, public Receiver<EntityDeactivatedEvent>, public Receiver<EntityActivatedEvent>
{
    struct DebugInfo
    {
        std::vector<VoxelRegion> lastRevoxelizationRegions;
    };
public:
    explicit VoxelizationPass();

    void init(float extentWorldLevel0, float virtualExtentWorldLevel0);

    void update() override;

    const DebugInfo& getDebugInfo() const { return m_debugInfo; }
private:
    void computeRevoxelizationRegionsClipmap(uint32_t clipmapLevel, std::vector<VoxelRegion> &clipRegions,
                                             std::vector<VoxelRegion> *revoxelizationRegions, const BBox &curBBox);

    /**
     * Compute the camera movement change delta (camera last pos and given clipmap region pos)
     * in voxel coordinates.
     */
    glm::ivec3 computeChangeDeltaV(uint32_t clipmapLevel, const VoxelRegion& clipRegions, const BBox& cameraRegionBBox);

    void computeRevoxelizationRegionsDynamicEntities(bool isVirtual, std::vector<VoxelRegion> &clipRegions, std::vector<VoxelRegion> *revoxelizationRegions);

    void receive(const EntityDeactivatedEvent& e) override;
    void receive(const EntityActivatedEvent& e) override;

    void recordDebugInfo();
private:
    std::shared_ptr<Shader> m_voxelizeShader;

    std::vector<VoxelRegion> m_revoxelizationRegions[CLIP_REGION_COUNT];
    std::vector<Entity> m_activatedDeactivatedEntities;
    std::vector<BBox> m_portionsOfDynamicEntities;
    std::vector<VoxelRegion> m_clipRegions;

    // A portion consisting of a multiple of voxel size is revoxelized
    int m_minChange[CLIP_REGION_COUNT] = {1, 1, 1, 1, 1, 1};

    Texture3D* m_voxelOpacity{nullptr};

    DebugInfo m_debugInfo;
    int m_overestimationWidth{ 0 };

    bool m_forceFullRevoxelization{ false };

    // Virtual object related variables
    Texture3D* m_virtualVoxelOpacity{nullptr};
    std::vector<VoxelRegion> m_virtualRevoxelizationRegions[VIRTUAL_CLIP_REGION_COUNT];
    std::vector<VoxelRegion> m_virtualClipRegions;
    bool m_isVirtualVoxelized{ false };
};
