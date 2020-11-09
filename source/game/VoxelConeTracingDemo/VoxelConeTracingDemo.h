#pragma once
#include <engine/Game.h>
#include <engine/input/Input.h>
#include <engine/resource/Model.h>
#include <engine/rendering/Material.h>
#include <engine/rendering/Texture3D.h>
#include <engine/util/Timer.h>
#include <engine/geometry/BBox.h>
#include "gui/VoxelConeTracingGUI.h"
#include "engine/rendering/voxelConeTracing/ClipmapUpdatePolicy.h"

class VoxelConeTracingDemo : public Game, InputHandler
{
public:
    VoxelConeTracingDemo();

    void update() override;
    void initUpdate() override;

    void moveCamera(Seconds deltaTime) const;

protected:
    void onKeyDown(SDL_Keycode keyCode) override;

private:
    void init3DVoxelTextures();

    BBox getBBox(size_t clipmapLevel, glm::vec3 center, float clipRegionBBoxExtentL0) const;

    void createDemoScene();
    void animateDirLight();
    void animateSphereRoughness();
    void animateCameraTransform();

    void updateCameraClipRegions();
    void updateVirtualClipRegions();

    ClipmapUpdatePolicy::Type getSelectedClipmapUpdatePolicyType() const;

private:
    std::unique_ptr<RenderPipeline> m_renderPipeline;
    std::vector<BBox> m_clipRegionBBoxes;
    std::vector<BBox> m_virtualClipRegionBBoxes;
    std::unique_ptr<ClipmapUpdatePolicy> m_clipmapUpdatePolicy;
    std::unique_ptr<ClipmapUpdatePolicy> m_virtualClipmapUpdatePolicy;

    // ClipRegion extent at level 0 - next level covers twice as much space as the previous level
    float m_clipRegionBBoxExtentL0{16.0f};
    float m_virtualClipRegionBBoxExtentL0{.0f};

    Texture3D m_voxelOpacity;
    Texture3D m_voxelRadiance;
    Texture3D m_voxelNormal;
    Texture3D m_voxelReflectance;

    Texture3D m_virtualVoxelOpacity;
    Texture3D m_virtualVoxelRadiance;
    Texture3D m_virtualVoxelNormal;
    Texture3D m_virtualVoxelDiffuse;
    Texture3D m_virtualVoxelSpecularA;
    glm::vec3 virtualMin, virtualMax;

    std::unique_ptr<VoxelConeTracingGUI> m_gui;
    bool m_guiEnabled{false};

    glm::vec3 m_scenePosition;
    Entity m_directionalLight;
    Entity m_sphere;
    ComponentPtr<Transform> virtualTransform;

    std::vector<glm::vec3> translations;
    std::vector<glm::quat> rotations;
    std::vector<Texture2D> backgroundImages;

// #define RISE103
// #define DASAN106
#define DASAN613

    const std::string cameraFilename = "../camera.txt";
#ifdef RISE103
    const std::string learningSceneDir = "../../neon/asset/rise103_learning6/";
    const std::string renderingSceneDir = "../../neon/asset/rise103_learning6/";
    const std::string sceneObjectFilename = "rise103_centered.obj";
    const std::string scenePCFilename = "cloud_rise103_learning_subsample=0.25cm_centered.ply";
    const std::string poseFilename = "../pose_rise103_learning6_centered.txt";
#elif defined(DASAN613)
    const std::string learningSceneDir = "../../neon/asset/dasan613_learning1/";
    const std::string renderingSceneDir = "../../neon/asset/dasan613_learning1/";
    // const std::string renderingSceneDir = "../../neon/asset/dasan613_rendering5-2/";
    const std::string sceneObjectFilename = "dasan613_tsdf3_centered.obj";
    const std::string scenePCFilename = "cloud_dasan613_learning_subsample=0.25cm_centered.ply";
    // const std::string poseFilename = "../pose_dasan613_rendering5-2_centered.txt";
    const std::string poseFilename = "../pose_dasan613_learning_centered.txt";
#elif defined(DASAN106)
    const std::string learningSceneDir = "../../neon/asset/dasan106_learning3/";
    const std::string renderingSceneDir = "../../neon/asset/dasan106_learning3/";
    // const std::string renderingSceneDir = "../../neon/asset/dasan106_rendering1/";
    const std::string sceneObjectFilename = "dasan106_centered.obj";
    const std::string scenePCFilename = "cloud_dasan106_learning_subsampled=0.25cm_centered.ply";
    // const std::string poseFilename = "../pose_dasan106_rendering1_centered.txt";
    const std::string poseFilename = "../pose_dasan106_learning3_centered.txt";
#endif
    const glm::vec3 centeringDasan106 = glm::vec3(4.407966, 5.205061, -4.503145);
    const glm::vec3 centeringDasan613 = glm::vec3(5.43885, 5.42074, - 5.40352);
    const glm::vec3 centeringRise103 = glm::vec3(5.015924, 4.872949, -5.523078);
};
