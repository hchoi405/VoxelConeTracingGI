#include "VoxelConeTracingDemo.h"
#include <engine/util/Random.h>
#include <engine/util/Timer.h>
#include <engine/rendering/util/GLUtil.h>
#include <engine/Engine.h>
#include <engine/rendering/Screen.h>
#include <engine/util/morton/morton.h>
#include <engine/util/util.h>
#include <engine/resource/ResourceManager.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <engine/rendering/architecture/RenderPipeline.h>
#include <engine/ecs/ECS.h>
#include <engine/rendering/lights/DirectionalLight.h>
#include "engine/rendering/voxelConeTracing/VoxelizationPass.h"
#include "engine/rendering/debug/DebugRenderer.h"
#include "engine/util/ECSUtil/ECSUtil.h"
#include "engine/rendering/renderPasses/SceneGeometryPass.h"
#include "engine/rendering/renderPasses/ShadowMapPass.h"
#include "engine/rendering/voxelConeTracing/RadianceInjectionPass.h"
#include "engine/rendering/voxelConeTracing/WrapBorderPass.h"
#include "engine/rendering/voxelConeTracing/GIPass.h"
#include "engine/rendering/voxelConeTracing/SphericalImagePass.h"
#include "engine/rendering/renderPasses/ForwardScenePass.h"
#include "engine/rendering/voxelConeTracing/settings/VoxelConeTracingSettings.h"
#include "engine/util/commands/RotationCommand.h"
#include "engine/util/commands/MaterialCommand.h"
#include "engine/util/commands/TransformCommand.h"
#include "engine/util/commands/CommandChain.h"
#include "engine/util/QueryManager.h"
#include "engine/rendering/renderer/MeshRenderers.h"
#include "engine/util/ECSUtil/EntityCreator.h"

using namespace glm;
#include "engine/rendering/tvcg17/common/commonStruct.h"

VoxelConeTracingDemo::VoxelConeTracingDemo() {
    Input::subscribe(this);

    Random::randomize();

    std::cout << "reading seuquence started... ";

    // Load pose
    std::ifstream in(poseFilename);
    if (!in.is_open()) {
        std::cout << "Failed to open " << poseFilename << std::endl;
        exit(-1);
    }
    while (!in.eof()) {
        glm::vec3 pos;
        glm::quat rot;
        in >> pos[0];
        in >> pos[1];
        in >> pos[2];
        in >> rot.x;
        in >> rot.y;
        in >> rot.z;
        in >> rot.w;
        translations.push_back(pos);
        rotations.push_back(rot);
    }
    in.close();
    std::cout << "finished!" << std::endl;
}

void VoxelConeTracingDemo::initUpdate() {
    ResourceManager::setShaderIncludePath("shaders");
    ResourceManager::setShaderIncludePath("../source/engine/rendering/tvcg17/common");

    createDemoScene();

    DebugRenderer::init();

    init3DVoxelTextures();

    // Load background
    int numSequence = translations.size();
    std::cout << "numSequence: " << numSequence << std::endl;
    if (numSequence > 0) {
        backgroundImages.resize(numSequence);
        for (int i = DEMO_SETTINGS.animateFrame.min; i <= DEMO_SETTINGS.animateFrame.max /* numSequence */; ++i) {
            std::stringstream ss;
            ss << renderingSceneDir << "color/";
            ss << std::setfill('0') << std::setw(5) << i;
            ss << ".png";
            std::cout << "load: " << ss.str() << std::endl;
            backgroundImages[i].load(ss.str());
        }
    }

    static Entity camera = ECS::getEntityByName("Camera");
    static auto camTransform = camera.getComponent<Transform>();
    camTransform->setPosition(translations[DEMO_SETTINGS.animateFrame]);
    camTransform->setRotation(rotations[DEMO_SETTINGS.animateFrame]);

    m_renderPipeline = std::make_unique<RenderPipeline>(MainCamera);
    m_gui = std::make_unique<VoxelConeTracingGUI>(m_renderPipeline.get());
    m_clipmapUpdatePolicy =
        std::make_unique<ClipmapUpdatePolicy>(ClipmapUpdatePolicy::Type::ONE_PER_FRAME_PRIORITY, CLIP_REGION_COUNT);
    m_virtualClipmapUpdatePolicy = std::make_unique<ClipmapUpdatePolicy>(
        ClipmapUpdatePolicy::Type::ONE_PER_FRAME_PRIORITY, VIRTUAL_CLIP_REGION_COUNT);
    auto sceneEntity = ECS::getEntityByName("name_scene");
    m_clipRegionBBoxExtentL0 = sceneEntity.getComponent<Transform>()->getBBox().maxExtent();
    std::cout << "Scene bbox: " << sceneEntity.getComponent<Transform>()->getBBox() << std::endl;
    std::cout << "m_clipRegionBBoxExtentL0: " << m_clipRegionBBoxExtentL0 << std::endl;
    m_virtualClipRegionBBoxExtentL0 = virtualTransform->getBBox().maxExtent() * 1.1f;
    std::cout << "m_virtualClipRegionBBoxExtentL0: " << m_virtualClipRegionBBoxExtentL0 << std::endl;
    virtualMin = virtualTransform->getBBox().min();
    virtualMax = virtualTransform->getBBox().max();
    m_renderPipeline->putPtr("virtualMin", &virtualMin);
    m_renderPipeline->putPtr("virtualMax", &virtualMax);

    // Set render pipeline input
    m_renderPipeline->putPtr("VoxelOpacity", &m_voxelOpacity);
    m_renderPipeline->putPtr("VoxelRadiance", &m_voxelRadiance);
    m_renderPipeline->putPtr("VoxelNormal", &m_voxelNormal);
    m_renderPipeline->putPtr("VoxelReflectance", &m_voxelReflectance);
    m_renderPipeline->putPtr("ClipRegionBBoxes", &m_clipRegionBBoxes);
    m_renderPipeline->putPtr("ClipmapUpdatePolicy", m_clipmapUpdatePolicy.get());

    m_renderPipeline->putPtr("VirtualVoxelOpacity", &m_virtualVoxelOpacity);
    m_renderPipeline->putPtr("VirtualVoxelRadiance", &m_virtualVoxelRadiance);
    m_renderPipeline->putPtr("VirtualVoxelNormal", &m_virtualVoxelNormal);
    m_renderPipeline->putPtr("VirtualVoxelDiffuse", &m_virtualVoxelDiffuse);
    m_renderPipeline->putPtr("VirtualVoxelSpecularA", &m_virtualVoxelSpecularA);
    m_renderPipeline->putPtr("VirtualClipRegionBBoxes", &m_virtualClipRegionBBoxes);
    m_renderPipeline->putPtr("VirtualClipmapUpdatePolicy", m_virtualClipmapUpdatePolicy.get());

    updateCameraClipRegions();
    updateVirtualClipRegions();

    // Add render passes to the pipeline
    m_renderPipeline->addRenderPasses(
        std::make_shared<SceneGeometryPass>(),
        std::make_shared<VoxelizationPass>()
        // ,std::make_shared<ShadowMapPass>(SHADOW_SETTINGS.shadowMapResolution) // Don't use the shadow map
        ,
        std::make_shared<RadianceInjectionPass>(), std::make_shared<WrapBorderPass>(), std::make_shared<GIPass>()
        // ,std::make_shared<SphericalImagePass>()
        //,std::make_shared<ForwardScenePass>()
    );

    // RenderPass initializations
    m_renderPipeline->getRenderPass<VoxelizationPass>()->init(m_clipRegionBBoxExtentL0,
                                                              m_virtualClipRegionBBoxExtentL0);
    std::cout << "m_clipRegionBBoxExtentL0: " << m_clipRegionBBoxExtentL0 << std::endl;

    // Deactivate after construction of VoxelizationPass to receive deactivated event
    // virtualTransform->getOwner().setActive(false);
    m_initializing = false;
}

void VoxelConeTracingDemo::update() {
    static bool once = true;

    // m_gui->selectEntity(virtualTransform->getOwner(), false);
    // if (once)
    //     m_gui->showVoxelVisualizationOptions();

    m_renderPipeline->put<GLuint>("BackgroundTexture", backgroundImages[DEMO_SETTINGS.animateFrame]);
    static Entity camera = ECS::getEntityByName("Camera");
    static auto camTransform = camera.getComponent<Transform>();
    if (DEMO_SETTINGS.animateCamera) {
        camTransform->setPosition(translations[DEMO_SETTINGS.animateFrame]);
        camTransform->setRotation(rotations[DEMO_SETTINGS.animateFrame]);
    }

    m_clipmapUpdatePolicy->setType(getSelectedClipmapUpdatePolicyType());
    m_clipmapUpdatePolicy->update();
    m_virtualClipmapUpdatePolicy->setType(getSelectedClipmapUpdatePolicyType());
    m_virtualClipmapUpdatePolicy->update();

    moveCamera(Time::deltaTime());

    glFrontFace(GL_CW);

    updateCameraClipRegions();
    updateVirtualClipRegions();

    for (auto c : ECS::getEntitiesWithComponents<CameraComponent>()) {
        c.getComponent<CameraComponent>()->updateViewMatrix();
    }

    bool giPipeline = RENDERING_SETTINGS.pipeline.asString() == "GI";
    bool forwardPipeline = RENDERING_SETTINGS.pipeline.asString() == "Forward";

    if (giPipeline) {
        // m_renderPipeline->getRenderPass<ForwardScenePass>()->setEnabled(false);
        m_renderPipeline->getRenderPass<GIPass>()->setEnabled(true);
        m_renderPipeline->getRenderPass<VoxelizationPass>()->setEnabled(once);
        // Use injection pass only once because one update is enough for point cloud
        m_renderPipeline->getRenderPass<RadianceInjectionPass>()->setEnabled(once);
        once = false;
        m_renderPipeline->getRenderPass<SceneGeometryPass>()->setEnabled(true);
    }

    if (forwardPipeline) {
        // m_renderPipeline->getRenderPass<ForwardScenePass>()->setEnabled(true);
        m_renderPipeline->getRenderPass<GIPass>()->setEnabled(false);
        m_renderPipeline->getRenderPass<VoxelizationPass>()->setEnabled(true);
        m_renderPipeline->getRenderPass<RadianceInjectionPass>()->setEnabled(false);
        m_renderPipeline->getRenderPass<SceneGeometryPass>()->setEnabled(false);
    }

    animateDirLight();
    animateSphereRoughness();
    animateCameraTransform();

    GL::setViewport(MainCamera->getViewport());

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    DebugRenderInfo info(MainCamera->view(), MainCamera->proj(), MainCamera->getPosition());
    DebugRenderer::begin(info);
    DebugRenderer::update();

    glEnable(GL_DEPTH_TEST);

    m_renderPipeline->update();

    ECS::lateUpdate();

    if (m_guiEnabled) m_gui->update();

    m_gui->onVoxelVisualization();

    DebugRenderer::end();
}

void VoxelConeTracingDemo::moveCamera(Seconds deltaTime) const {
    float speed = DEMO_SETTINGS.cameraSpeed * deltaTime;

    if (Input::isKeyDown(SDL_SCANCODE_W)) MainCamera->walk(speed);

    if (Input::isKeyDown(SDL_SCANCODE_S)) MainCamera->walk(-speed);

    if (Input::isKeyDown(SDL_SCANCODE_A)) MainCamera->strafe(-speed);

    if (Input::isKeyDown(SDL_SCANCODE_D)) MainCamera->strafe(speed);

    if (Input::isKeyDown(SDL_SCANCODE_Q)) MainCamera->elevate(-speed);

    if (Input::isKeyDown(SDL_SCANCODE_E)) MainCamera->elevate(speed);

    if (Input::rightDrag().isDragging()) {
        float dx = -Input::rightDrag().getDragDelta().x;
        float dy = Input::rightDrag().getDragDelta().y;

        MainCamera->pitch(math::toRadians(dy * 0.1f));
        MainCamera->rotateY(math::toRadians(dx * 0.1f));
        // SDL_SetRelativeMouseMode(SDL_TRUE);
    } else {
        SDL_SetRelativeMouseMode(SDL_FALSE);
    }
}

void VoxelConeTracingDemo::onKeyDown(SDL_Keycode keyCode) {
    switch (keyCode) {
        case SDLK_F1:
            m_guiEnabled = !m_guiEnabled;
            break;
        case SDLK_c: {
            // TODO: check where this error comes from
            GL_ERROR_CHECK();
            auto sgPass = m_renderPipeline->getRenderPass<SceneGeometryPass>();
            auto texture = sgPass->getRenderTexturePtr(GL_COLOR_ATTACHMENT5);
            texture->save("test", true, sgPass->getFBO());
            break;
        }
        case SDLK_f: {
            auto camTransform = MainCamera->getComponent<Transform>();
            std::ofstream cameraFile(cameraFilename);
            auto euler = glm::degrees(camTransform->getEulerAngles());
            cameraFile << camTransform->getPosition()[0] << " " << camTransform->getPosition()[1] << " "
                       << camTransform->getPosition()[2] << std::endl
                       << euler[0] << " " << euler[1] << " " << euler[2];
            cameraFile.close();
            std::cout << "Saved camera position " << camTransform->getPosition() << " and rotation " << euler << " to "
                      << cameraFilename << std::endl;
            break;
        }
        case SDLK_r:
            MainCamera->rotateY(math::toRadians(180));
            break;
        case SDLK_z:
            m_gui->selectEntity(virtualTransform->getOwner(), false);
            break;
        case SDLK_F5:
            m_engine->requestScreenshot();
            break;
        default:
            break;
    }
}

void VoxelConeTracingDemo::init3DVoxelTextures() {
    // GLint filter = GL_LINEAR;
    GLint filter = GL_NEAREST;
    GLint wrapS = GL_CLAMP_TO_BORDER;
    GLint wrapT = GL_CLAMP_TO_BORDER;
    GLint wrapR = GL_CLAMP_TO_BORDER;

    // Each anisotropic voxel needs multiple values (per face = 6). Furthermore multiple clipmap regions are required.
    // Per attribute everything is stored in one 3D texture:
    // In X Direction -> voxelFace0,voxelFace1,...,voxelFace5
    // In Y Direction -> clipmap L0,L1,...,Ln
    // So the position in [0, resolution - 1]^3 of a specific inner texture is accessed by:
    // samplePos(x, y, z, faceIndex, clipmapLevel) = (x, y, z) + (faceIndex, clipmapLevel, 0) * resolution
    // To ensure correct interpolation at borders we thus need to extend the resolution of each inner texture
    // by 2 in each dimension and copy in a dedicated render pass (WrapBorderPass) the border of the other side
    // to get GL_REPEAT as the texture wrapping mode.
    GLsizei resolutionWithBorder = VOXEL_RESOLUTION + 2;

    m_voxelOpacity.create(resolutionWithBorder * FACE_COUNT, CLIP_REGION_COUNT * resolutionWithBorder,
                          resolutionWithBorder, GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE, Texture3DSettings::Custom);
    m_voxelOpacity.bind();
    m_voxelOpacity.setParameteri(GL_TEXTURE_WRAP_S, wrapS);
    m_voxelOpacity.setParameteri(GL_TEXTURE_WRAP_T, wrapT);
    m_voxelOpacity.setParameteri(GL_TEXTURE_WRAP_R, wrapR);
    m_voxelOpacity.setParameteri(GL_TEXTURE_MIN_FILTER, filter);
    m_voxelOpacity.setParameteri(GL_TEXTURE_MAG_FILTER, filter);

    m_voxelRadiance.create(resolutionWithBorder * FACE_COUNT, CLIP_REGION_COUNT * resolutionWithBorder,
                           resolutionWithBorder, GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE, Texture3DSettings::Custom);
    m_voxelRadiance.bind();
    m_voxelRadiance.setParameteri(GL_TEXTURE_WRAP_S, wrapS);
    m_voxelRadiance.setParameteri(GL_TEXTURE_WRAP_T, wrapT);
    m_voxelRadiance.setParameteri(GL_TEXTURE_WRAP_R, wrapR);
    m_voxelRadiance.setParameteri(GL_TEXTURE_MIN_FILTER, filter);
    m_voxelRadiance.setParameteri(GL_TEXTURE_MAG_FILTER, filter);

    m_voxelNormal.create(resolutionWithBorder * FACE_COUNT, CLIP_REGION_COUNT * resolutionWithBorder,
                         resolutionWithBorder, GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE, Texture3DSettings::Custom);
    m_voxelNormal.bind();
    m_voxelNormal.setParameteri(GL_TEXTURE_WRAP_S, wrapS);
    m_voxelNormal.setParameteri(GL_TEXTURE_WRAP_T, wrapT);
    m_voxelNormal.setParameteri(GL_TEXTURE_WRAP_R, wrapR);
    m_voxelNormal.setParameteri(GL_TEXTURE_MIN_FILTER, filter);
    m_voxelNormal.setParameteri(GL_TEXTURE_MAG_FILTER, filter);

    m_voxelReflectance.create(resolutionWithBorder * FACE_COUNT, CLIP_REGION_COUNT * resolutionWithBorder,
                              resolutionWithBorder, GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE, Texture3DSettings::Custom);
    m_voxelReflectance.bind();
    m_voxelReflectance.setParameteri(GL_TEXTURE_WRAP_S, wrapS);
    m_voxelReflectance.setParameteri(GL_TEXTURE_WRAP_T, wrapT);
    m_voxelReflectance.setParameteri(GL_TEXTURE_WRAP_R, wrapR);
    m_voxelReflectance.setParameteri(GL_TEXTURE_MIN_FILTER, filter);
    m_voxelReflectance.setParameteri(GL_TEXTURE_MAG_FILTER, filter);

    resolutionWithBorder = VIRTUAL_VOXEL_RESOLUTION + 2;

    m_virtualVoxelOpacity.create(resolutionWithBorder * FACE_COUNT, VIRTUAL_CLIP_REGION_COUNT * resolutionWithBorder,
                                 resolutionWithBorder, GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE, Texture3DSettings::Custom);
    m_virtualVoxelOpacity.bind();
    m_virtualVoxelOpacity.setParameteri(GL_TEXTURE_WRAP_S, wrapS);
    m_virtualVoxelOpacity.setParameteri(GL_TEXTURE_WRAP_T, wrapT);
    m_virtualVoxelOpacity.setParameteri(GL_TEXTURE_WRAP_R, wrapR);
    m_virtualVoxelOpacity.setParameteri(GL_TEXTURE_MIN_FILTER, filter);
    m_virtualVoxelOpacity.setParameteri(GL_TEXTURE_MAG_FILTER, filter);

    m_virtualVoxelRadiance.create(resolutionWithBorder * FACE_COUNT, VIRTUAL_CLIP_REGION_COUNT * resolutionWithBorder,
                                  resolutionWithBorder, GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE, Texture3DSettings::Custom);
    m_virtualVoxelRadiance.bind();
    m_virtualVoxelRadiance.setParameteri(GL_TEXTURE_WRAP_S, wrapS);
    m_virtualVoxelRadiance.setParameteri(GL_TEXTURE_WRAP_T, wrapT);
    m_virtualVoxelRadiance.setParameteri(GL_TEXTURE_WRAP_R, wrapR);
    m_virtualVoxelRadiance.setParameteri(GL_TEXTURE_MIN_FILTER, filter);
    m_virtualVoxelRadiance.setParameteri(GL_TEXTURE_MAG_FILTER, filter);

    m_virtualVoxelDiffuse.create(resolutionWithBorder * FACE_COUNT, VIRTUAL_CLIP_REGION_COUNT * resolutionWithBorder,
                                 resolutionWithBorder, GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE, Texture3DSettings::Custom);
    m_virtualVoxelDiffuse.bind();
    m_virtualVoxelDiffuse.setParameteri(GL_TEXTURE_WRAP_S, wrapS);
    m_virtualVoxelDiffuse.setParameteri(GL_TEXTURE_WRAP_T, wrapT);
    m_virtualVoxelDiffuse.setParameteri(GL_TEXTURE_WRAP_R, wrapR);
    m_virtualVoxelDiffuse.setParameteri(GL_TEXTURE_MIN_FILTER, filter);
    m_virtualVoxelDiffuse.setParameteri(GL_TEXTURE_MAG_FILTER, filter);

    m_virtualVoxelSpecularA.create(resolutionWithBorder * FACE_COUNT, VIRTUAL_CLIP_REGION_COUNT * resolutionWithBorder,
                                   resolutionWithBorder, GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE,
                                   Texture3DSettings::Custom);
    m_virtualVoxelSpecularA.bind();
    m_virtualVoxelSpecularA.setParameteri(GL_TEXTURE_WRAP_S, wrapS);
    m_virtualVoxelSpecularA.setParameteri(GL_TEXTURE_WRAP_T, wrapT);
    m_virtualVoxelSpecularA.setParameteri(GL_TEXTURE_WRAP_R, wrapR);
    m_virtualVoxelSpecularA.setParameteri(GL_TEXTURE_MIN_FILTER, filter);
    m_virtualVoxelSpecularA.setParameteri(GL_TEXTURE_MAG_FILTER, filter);

    m_virtualVoxelNormal.create(resolutionWithBorder * FACE_COUNT, VIRTUAL_CLIP_REGION_COUNT * resolutionWithBorder,
                                resolutionWithBorder, GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE, Texture3DSettings::Custom);
    m_virtualVoxelNormal.bind();
    m_virtualVoxelNormal.setParameteri(GL_TEXTURE_WRAP_S, wrapS);
    m_virtualVoxelNormal.setParameteri(GL_TEXTURE_WRAP_T, wrapT);
    m_virtualVoxelNormal.setParameteri(GL_TEXTURE_WRAP_R, wrapR);
    m_virtualVoxelNormal.setParameteri(GL_TEXTURE_MIN_FILTER, filter);
    m_virtualVoxelNormal.setParameteri(GL_TEXTURE_MAG_FILTER, filter);

    GL_ERROR_CHECK();
}

/**
 * return: bounding box of clipmapLevel in world coordinate around center
 */
BBox VoxelConeTracingDemo::getBBox(size_t clipmapLevel, glm::vec3 center, float clipRegionBBoxExtentL0) const {
    float halfSize = 0.5f * clipRegionBBoxExtentL0 * std::exp2f(float(clipmapLevel));

    return BBox(center - halfSize, center + halfSize);
}

void VoxelConeTracingDemo::createDemoScene() {
    m_scenePosition = glm::vec3(0.0f);

    Entity camera = ECS::createEntity("Camera");
    camera.addComponent<Transform>();
    camera.addComponent<CameraComponent>();

    auto camComponent = camera.getComponent<CameraComponent>();
    auto camTransform = camera.getComponent<Transform>();

    MainCamera = camComponent;

    camComponent->setPerspective(48.8694, float(Screen::getWidth()), float(Screen::getHeight()), 0.3f, 30.0f);

    std::ifstream cameraFile(cameraFilename);
    glm::vec3 cameraPositionOffset;
    glm::vec3 cameraRotation;
    if (!cameraFile.is_open()) {
        std::cout << "failed to read camera.txt" << std::endl;
    } else {
        if (rotations.empty()) {
            std::cout << "Load camera.txt" << std::endl;
            cameraFile >> cameraPositionOffset[0] >> cameraPositionOffset[1] >> cameraPositionOffset[2];
            cameraFile >> cameraRotation[0] >> cameraRotation[1] >> cameraRotation[2];
            camTransform->setEulerAngles(glm::radians(cameraRotation));
            camTransform->setPosition(m_scenePosition + cameraPositionOffset);
            // Set rotation manually
            camTransform->lookAt(vec3(0.85f, 0.5f, -0.5f));
        } else {
            std::cout << "Set camera's parameter using the first element of sequence" << std::endl;
            camTransform->setPosition(translations[0]);
            camTransform->setRotation(rotations[0]);
        }
    }

    m_engine->registerCamera(camComponent);

    auto shader = ResourceManager::getShader("shaders/forwardShadingPass.vert", "shaders/forwardShadingPass.frag",
                                             {"in_pos", "in_normal", "in_tangent", "in_bitangent", "in_uv"});

    ResourceManager::getModel(learningSceneDir + sceneObjectFilename)->name = "name_scene";
    auto sceneRootEntity = ECSUtil::loadMeshEntities(learningSceneDir + sceneObjectFilename, shader, learningSceneDir,
                                                     glm::vec3(1.f), false);
    // sceneRootEntity->setEulerAngles(glm::vec3(math::toRadians(90.f), math::toRadians(0.f), math::toRadians(0.f)));
    std::cout << "min: " << sceneRootEntity->getBBox().min() << std::endl;
    std::cout << "max: " << sceneRootEntity->getBBox().max() << std::endl;

    // Point Clout Entity
    std::string pcEntityName = "PointCloud";
    std::string cloudFilename = learningSceneDir + scenePCFilename;

    // Load point cloud before loadMeshEntities to give a name
    // This is ok because loadMeshEntities() first tries to find
    // whether the model given with the filename is already loaded
    ResourceManager::getModel(cloudFilename)->name = pcEntityName;
    auto pcEntityTransform = ECSUtil::loadMeshEntities(cloudFilename, shader, learningSceneDir, glm::vec3(1.f), false);
    pcEntityTransform->getComponent<MeshRenderer>()->getMesh()->setRenderMode(GL_POINTS, 0);
    // pcEntityTransform->setEulerAngles(glm::vec3(math::toRadians(90.f), math::toRadians(0.f), math::toRadians(0.f)));
    pcEntityTransform->setPosition(glm::vec3(m_scenePosition));

    auto pcEntity = ECS::getEntityByName(pcEntityName);
    pcEntity.setActive(false);

    // Virtual sphere
    // m_sphere = EntityCreator::createSphere("virtualObject", glm::vec3(0), glm::vec3(1.f));
    std::string virtualObjectDir = "../../neon/asset/mesh/";
    std::string virtualEntityName = "virtualObject";

    {
#ifdef DASAN613
    /*   // Lucy
      std::string virtualObjectFilename = "lucy_co_tri_simplified4_centered.ply";
      ResourceManager::getModel(virtualObjectDir + virtualObjectFilename)->name = virtualEntityName;
      virtualTransform = ECSUtil::loadMeshEntities(virtualObjectDir + virtualObjectFilename, shader, virtualObjectDir,
                                                   glm::vec3(1.f), false);
      virtualTransform->setPosition(
          // Original object center
          glm::vec3(5.23f, 5.90143f, -4.57f)
           - centeringDasan613
      );  // Dasan613

      auto sphereMaterial = EntityCreator::createMaterial();
      sphereMaterial->setColor("u_color", glm::vec4(0.48f, 0.392f, 0.114, 1.0f));
      virtualTransform->getOwner().getComponent<MeshRenderer>()->setMaterial(sphereMaterial, 0);
      // Virtual should be set at the root entity
      m_sphere = ECS::getEntityByName(virtualEntityName);
      m_sphere.setVirtual(true);
      m_sphere.setActive(true); */
#endif
    }

    {
        /*   // Cube
          std::string virtualObjectFilename = "cube_centered.obj";
          ResourceManager::getModel(virtualObjectDir + virtualObjectFilename)->name = virtualEntityName;
          virtualTransform = ECSUtil::loadMeshEntities(virtualObjectDir + virtualObjectFilename, shader,
          virtualObjectDir,
                                                       glm::vec3(1.f), true);
          // virtualTransform->setPosition(glm::vec3(5.4f, 6.24f, -4.57f)); // Dasan613
          virtualTransform->setPosition(glm::vec3(6.2384f, 6.04f, -5.33125f) - centeringDasan106 + glm::vec3(0, 0.25,
          0)); virtualTransform->setLocalScale(glm::vec3(0.5f));

          // Plate
          std::string virtualObjectFilename = "plate_centered.obj";
          ResourceManager::getModel(virtualObjectDir + virtualObjectFilename)->name = virtualEntityName;
          virtualTransform = ECSUtil::loadMeshEntities(virtualObjectDir + virtualObjectFilename, shader,
          virtualObjectDir,
                                                       glm::vec3(1.f), true);
          virtualTransform->setPosition(glm::vec3(3.950f, 4.55f, -4.570) - centeringDasan106 +
                                        glm::vec3(0, 0.5, 0));
          auto rot = glm::rotate(glm::radians(-45.f), glm::vec3(0, 1, 0));
          virtualTransform->setRotation(glm::toQuat(rot));
   */
    }

    {
//         // Sphere
//         std::string virtualObjectFilename = "sphere_high_centered.obj";
//         ResourceManager::getModel(virtualObjectDir + virtualObjectFilename)->name = virtualEntityName;
//         virtualTransform = ECSUtil::loadMeshEntities(virtualObjectDir + virtualObjectFilename, shader, virtualObjectDir,
//                                                      glm::vec3(5.f), false);
// #ifdef DASAN106
//         virtualTransform->setPosition(glm::vec3(6.2384f, 6.16864f, -5.33125f) - centeringDasan106 +
//                                       glm::vec3(0, 0.0125, 0));
// #elif defined(DASAN613)
//         //  ne::vec3f(5.4f, 6.35f, 4.57f)
//         virtualTransform->setPosition(glm::vec3(5.4f, 6.35f, -4.57f) - centeringDasan613);
//         std::cout << virtualTransform->getPosition() << std::endl;
// #endif

//         // Virtual should be set at the root entity
//         m_sphere = ECS::getEntityByName(virtualEntityName);
//         m_sphere.setVirtual(true);
//         m_sphere.setActive(true);
    }

    {
        // Signle Buddha
        std::string virtualObjectFilename = "buddha_centered.obj";
        ResourceManager::getModel(virtualObjectDir + virtualObjectFilename)->name = virtualEntityName;
        virtualTransform = ECSUtil::loadMeshEntities(virtualObjectDir + virtualObjectFilename, shader, virtualObjectDir,
                                                     glm::vec3(1.f), false);
#ifdef RISE103
        auto calcPos = glm::vec3(5.56768, 5.55787, -3.98624) - centeringRise103;
        // - glm::vec3(1.38286, 1.32292, 2.10885);
        virtualTransform->setPosition(calcPos);  // Rise103
        std::cout << "calcPos: " << calcPos << std::endl;
#endif

        auto sphereMaterial = EntityCreator::createMaterial();
        sphereMaterial->setFloat("u_shininess", 255.0f);
        sphereMaterial->setColor("u_color", glm::vec4(0.880392f, 0.768627f, 0.323725f, 1.0f));
        sphereMaterial->setColor("u_emissionColor", glm::vec3(0.0f));
        sphereMaterial->setColor("u_specularColor", glm::vec3(1.f));
        auto childTransform = virtualTransform->getChildren()[0];
        childTransform->getComponent<MeshRenderer>()->setMaterial(sphereMaterial, 0);
        // Virtual should be set at the root entity
        m_sphere = ECS::getEntityByName(virtualEntityName);
        m_sphere.setVirtual(true);
        m_sphere.setActive(false);
    }

    {
     /*     // sphere, bunny, lucy
         auto vo1 = ResourceManager::getModel(virtualObjectDir + "seq4_sphere_glass_2_centered.ply");
         vo1->name = "vo1";
         auto vo2 = ResourceManager::getModel(virtualObjectDir + "seq4_bunny_diffuse_2_centered.ply");
         vo2->name = "vo2";
         auto vo3 = ResourceManager::getModel(virtualObjectDir + "seq4_lucy_glossy_2_centered.ply");
         vo3->name = "vo3";

         auto virtualObject = ECS::createEntity("virtualObject");
         virtualObject.addComponent<Transform>();
         auto parentTransform = virtualObject.getComponent<Transform>();

         auto voTransform1 = ECSUtil::loadMeshEntities(vo1.get(), shader, "", glm::vec3(1.f), false);
         auto sphereMaterial = EntityCreator::createMaterial();
         sphereMaterial->setColor("u_color", glm::vec4(0.5f));  // Material flag
         voTransform1->getOwner().getComponent<MeshRenderer>()->setMaterial(sphereMaterial, 0);
         voTransform1->getOwner().setVirtual(true);
         voTransform1->getOwner().setActive(true);
         voTransform1->setParent(parentTransform);

         auto voTransform2 = ECSUtil::loadMeshEntities(vo2.get(), shader, "", glm::vec3(1.f), false);
         auto bunnyMaterial = EntityCreator::createMaterial();
         bunnyMaterial->setColor("u_color", glm::vec4(0.9f, 0.2f, 0.113725f, 1));  // Material flag
         voTransform2->getOwner().getComponent<MeshRenderer>()->setMaterial(bunnyMaterial, 0);
         voTransform2->getOwner().setVirtual(true);
         voTransform2->getOwner().setActive(true);
         voTransform2->setParent(parentTransform);

         auto voTransform3 = ECSUtil::loadMeshEntities(vo3.get(), shader, "", glm::vec3(1.f), false);
         auto lucyMaterial = EntityCreator::createMaterial();
         lucyMaterial->setColor("u_color", glm::vec4(0.48f, 0.392f, 0.114f, 1));  // Material flag
         voTransform3->getOwner().getComponent<MeshRenderer>()->setMaterial(lucyMaterial, 0);
         voTransform3->getOwner().setVirtual(true);
         voTransform3->getOwner().setActive(true);
         voTransform3->setParent(parentTransform);

         BBox parentBBox = voTransform1->getBBox();
         parentBBox.unite(voTransform2->getBBox());
         parentBBox.unite(voTransform3->getBBox());
         parentTransform->setBBox(parentBBox);

         // auto calcPos = glm::vec3(5.37086533, 6.20498467, -4.756995) - centeringDasan613; // v1
         auto calcPos = glm::vec3(5.3362, 6.3235, -4.57801) - centeringDasan613; // v2

         parentTransform->setPosition(calcPos);
         std::cout << "calculated pos: " << calcPos << std::endl;
         std::cout << "parent bbox: " << parentTransform->getBBox() << std::endl;
         virtualTransform = parentTransform; */
    }

    {
      /*    // Virtual Buddha
         auto vo1 = ResourceManager::getModel(virtualObjectDir + "buddha1_dasan106_simplified.ply");
         vo1->name = "vo1";
         auto vo2 = ResourceManager::getModel(virtualObjectDir + "buddha2_dasan106_simplified.ply");
         vo2->name = "vo2";

         auto virtualObject = ECS::createEntity("virtualObject");
         virtualObject.addComponent<Transform>();
         auto parentTransform = virtualObject.getComponent<Transform>();

         // buddha 1
         auto voTransform1 = ECSUtil::loadMeshEntities(vo1.get(), shader, "", glm::vec3(1.f), true);
         std::cout << "1 bbox: " << voTransform1->getBBox() << std::endl;
         // auto voTransform1 = EntityCreator::createSphere("virtualObject", glm::vec3(0),
         // glm::vec3(1.f)).getComponent<Transform>();
         auto buddhaMaterial1 = EntityCreator::createMaterial();
         buddhaMaterial1->setFloat("u_shininess", 255.0);
         buddhaMaterial1->setColor("u_color", glm::vec4(1, 0, 0, 1));
         buddhaMaterial1->setColor("u_emissionColor", glm::vec3(0.0f));
         buddhaMaterial1->setColor("u_specularColor", glm::vec3(1.f));
         voTransform1->getOwner().getComponent<MeshRenderer>()->setMaterial(buddhaMaterial1, 0);
         voTransform1->getOwner().setVirtual(true);
         voTransform1->getOwner().setActive(true);
         voTransform1->setParent(parentTransform);

         // buddha 2
         auto voTransform2 = ECSUtil::loadMeshEntities(vo2.get(), shader, "", glm::vec3(1.f), true);
         std::cout << "2 bbox: " << voTransform2->getBBox() << std::endl;
         auto buddhaMaterial2 = EntityCreator::createMaterial();
         buddhaMaterial2->setFloat("u_shininess", 255.0);
         buddhaMaterial2->setColor("u_color", glm::vec4(1));
         buddhaMaterial2->setColor("u_emissionColor", glm::vec3(0.0f));
         buddhaMaterial2->setColor("u_specularColor", glm::vec3(1.f));
         voTransform2->getOwner().getComponent<MeshRenderer>()->setMaterial(buddhaMaterial2, 0);
         voTransform2->getOwner().setVirtual(true);
         voTransform2->getOwner().setActive(true);
         // voTransform2->setLocalPosition(glm::vec3(1, -0.5, 0));
         // voTransform2->setLocalEulerAngles(glm::radians(glm::vec3(180.f, 0.f, 0.f)));
         voTransform2->setParent(parentTransform);

         // parentTransform->setPosition(glm::vec3(5.4f, 7.05f, -4.57f));
         // parentTransform->setPosition();
         // parentTransform->setEulerAngles(glm::radians(vec3(0, 0, 180.f)));
         parentTransform->getOwner().setVirtual(true);

         BBox parentBBox = voTransform1->getBBox();
         parentBBox.unite(voTransform2->getBBox());
         parentTransform->setBBox(parentBBox);

 #ifdef DASAN106
         auto calcPos = glm::vec3(6.130078, 5.9573895, -5.05876) - centeringDasan106;
         // parentTransform->setPosition(glm::vec3(1.8f, 0.7f, -0.75f));
         parentTransform->setPosition(calcPos);
         std::cout << "calculated pos: " << calcPos << std::endl;
 #endif

         std::cout << "parent bbox: " << parentTransform->getBBox() << std::endl;

         virtualTransform = parentTransform; */
    }

    if (sceneRootEntity) sceneRootEntity->setPosition(glm::vec3(m_scenePosition));
}

void VoxelConeTracingDemo::animateDirLight() {
    if (m_directionalLight && DEMO_SETTINGS.animateLight) {
        auto lightTransform = m_directionalLight.getComponent<Transform>();
        glm::vec3 lightPos = lightTransform->getPosition();
        glm::quat startRotation = glm::quat(
            glm::transpose(glm::lookAt(lightPos, glm::vec3(0.0f, 0.0f, -10.0f), glm::vec3(0.0f, 1.0f, 0.0f))));
        glm::quat dstRotation =
            glm::quat(glm::transpose(glm::lookAt(lightPos, glm::vec3(0.0f, 0.0f, 10.0f), glm::vec3(0.0f, 1.0f, 0.0f))));

        static std::shared_ptr<RotationCommand> rotationCommand =
            std::make_shared<RotationCommand>(lightTransform, startRotation, dstRotation, 10.0f);
        static std::shared_ptr<RotationCommand> rotationCommand2 =
            std::make_shared<RotationCommand>(lightTransform, dstRotation, startRotation, 10.0f);
        static CommandChain commandChain({rotationCommand, rotationCommand2}, true);

        commandChain(Time::deltaTime());
    }
}

void VoxelConeTracingDemo::animateSphereRoughness() {
    if (m_sphere && DEMO_SETTINGS.animateSphere) {
        auto sphereTransform = m_sphere.getComponent<Transform>();
        static glm::quat initialRot = sphereTransform->getRotation();

        static float startY = 0.506;
        static float endY = 0.8;

        std::shared_ptr<TransformCommand> tCommand = std::make_shared<TransformCommand>(
            sphereTransform, glm::vec3(1.35, startY, -1.3), glm::vec3(1.35, endY, -1.3), initialRot, initialRot,
            DEMO_SETTINGS.cameraSpeed);
        std::shared_ptr<TransformCommand> tCommand2 = std::make_shared<TransformCommand>(
            sphereTransform, glm::vec3(1.35, endY, -1.3), glm::vec3(1.35, startY, -1.3), initialRot, initialRot,
            DEMO_SETTINGS.cameraSpeed);
        static CommandChain commandChain({tCommand, tCommand2}, true);

        commandChain(Time::deltaTime());
    }
}

void VoxelConeTracingDemo::animateCameraTransform() {
    // static size_t frame = 0;
    if (DEMO_SETTINGS.animateCamera) {
        auto cameraTransform = MainCamera.getOwner().getComponent<Transform>();
        static glm::vec3 initialPos = cameraTransform->getPosition();
        static glm::quat initialRot = cameraTransform->getRotation();

        static std::shared_ptr<TransformCommand> tCommand = std::make_shared<TransformCommand>(
            cameraTransform, translations[DEMO_SETTINGS.animateFrame], translations[DEMO_SETTINGS.animateFrame + 1],
            rotations[DEMO_SETTINGS.animateFrame], rotations[DEMO_SETTINGS.animateFrame + 1],
            0.0333f / DEMO_SETTINGS.cameraSpeed);
        static CommandChain commandChain({tCommand}, false);
        if (commandChain.done()) {
            if (size_t(DEMO_SETTINGS.animateFrame.value) < rotations.size() - 1) {
                m_engine->requestScreenshot();
                m_guiEnabled = false;
                DEMO_SETTINGS.animateFrame.value++;
                tCommand = std::make_shared<TransformCommand>(
                    cameraTransform, translations[DEMO_SETTINGS.animateFrame],
                    translations[DEMO_SETTINGS.animateFrame + 1], rotations[DEMO_SETTINGS.animateFrame],
                    rotations[DEMO_SETTINGS.animateFrame + 1], 0.0333f / DEMO_SETTINGS.cameraSpeed);
                commandChain = CommandChain({tCommand}, false);
            } else {
                std::cout << "Go to initial" << std::endl;
                // Set to start position and location
                cameraTransform->setPosition(initialPos);
                cameraTransform->setRotation(initialRot);
                DEMO_SETTINGS.animateCamera.value = false;
                exit(0);
            }
        }

        commandChain(Time::deltaTime());
    } else {
        // DEMO_SETTINGS.animateFrame.value = 0;
    }
}

void VoxelConeTracingDemo::updateCameraClipRegions() {
    m_clipRegionBBoxes.clear();
    auto sceneEntity = ECS::getEntityByName("name_scene");
    // glm::vec3 center = MainCamera->getPosition();
    glm::vec3 center = sceneEntity.getComponent<Transform>()->getPosition();

    // Align with virtual object's clip region
    // glm::vec3 center = virtualTransform->getPosition();

    for (size_t i = 0; i < CLIP_REGION_COUNT; ++i)
        m_clipRegionBBoxes.push_back(getBBox(i, center, m_clipRegionBBoxExtentL0));
}

void VoxelConeTracingDemo::updateVirtualClipRegions() {
    // add bbox for each clip level i around the virtual object
    m_virtualClipRegionBBoxes.clear();

    glm::vec3 center = virtualTransform->getPosition();

    for (size_t i = 0; i < VIRTUAL_CLIP_REGION_COUNT; ++i) {
        m_virtualClipRegionBBoxes.push_back(getBBox(i, center, m_virtualClipRegionBBoxExtentL0));
    }
}

ClipmapUpdatePolicy::Type VoxelConeTracingDemo::getSelectedClipmapUpdatePolicyType() const {
    if (GI_SETTINGS.updateOneClipLevelPerFrame) return ClipmapUpdatePolicy::Type::ONE_PER_FRAME_PRIORITY;

    return ClipmapUpdatePolicy::Type::ALL_PER_FRAME;
}
