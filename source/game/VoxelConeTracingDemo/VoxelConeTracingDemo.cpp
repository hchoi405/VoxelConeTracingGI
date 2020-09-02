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

    std::ifstream in("../output.txt");
    if (!in.is_open()) {
        std::cout << "Failed to open output.txt" << std::endl;
        exit(-1);
    }
    std::cout << "reading started... ";
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

    m_renderPipeline = std::make_unique<RenderPipeline>(MainCamera);
    m_gui = std::make_unique<VoxelConeTracingGUI>(m_renderPipeline.get());
    m_clipmapUpdatePolicy =
        std::make_unique<ClipmapUpdatePolicy>(ClipmapUpdatePolicy::Type::ONE_PER_FRAME_PRIORITY, CLIP_REGION_COUNT);
    m_virtualClipmapUpdatePolicy = std::make_unique<ClipmapUpdatePolicy>(
        ClipmapUpdatePolicy::Type::ONE_PER_FRAME_PRIORITY, VIRTUAL_CLIP_REGION_COUNT);
    auto sceneEntity = ECS::getEntityByName("name_scene");
    m_clipRegionBBoxExtentL0 = sceneEntity.getComponent<Transform>()->getBBox().maxExtent() * 3;
    std::cout << "Scene bbox: " << sceneEntity.getComponent<Transform>()->getBBox() << std::endl;
    std::cout << "m_clipRegionBBoxExtentL0: " << m_clipRegionBBoxExtentL0 << std::endl;
    m_virtualClipRegionBBoxExtentL0 = virtualTransform->getBBox().maxExtent() * 1.5;
    std::cout << "m_virtualClipRegionBBoxExtentL0: " << m_virtualClipRegionBBoxExtentL0 << std::endl;

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
        std::make_shared<VoxelizationPass>(),
        // std::make_shared<ShadowMapPass>(SHADOW_SETTINGS.shadowMapResolution), // Don't use the shadow map
        std::make_shared<RadianceInjectionPass>(),
        std::make_shared<WrapBorderPass>(),
        std::make_shared<GIPass>(),
        // std::make_shared<SphericalImagePass>(),
        std::make_shared<ForwardScenePass>());

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
        m_renderPipeline->getRenderPass<ForwardScenePass>()->setEnabled(false);
        m_renderPipeline->getRenderPass<GIPass>()->setEnabled(true);
        m_renderPipeline->getRenderPass<VoxelizationPass>()->setEnabled(once);
        // Use injection pass only once because one update is enough for point cloud
        m_renderPipeline->getRenderPass<RadianceInjectionPass>()->setEnabled(once);
        once = false;
        m_renderPipeline->getRenderPass<SceneGeometryPass>()->setEnabled(true);
    }

    if (forwardPipeline) {
        m_renderPipeline->getRenderPass<ForwardScenePass>()->setEnabled(true);
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

    camTransform->setPosition(translations[340]);
    camTransform->setRotation(rotations[340]);
    m_engine->registerCamera(camComponent);

    auto shader = ResourceManager::getShader("shaders/forwardShadingPass.vert", "shaders/forwardShadingPass.frag",
                                             {"in_pos", "in_normal", "in_tangent", "in_bitangent", "in_uv"});

    ResourceManager::getModel("../../neon/asset/dasan613_learning1/dasan613_tsdf3.obj")->name = "name_scene";
    auto sceneRootEntity = ECSUtil::loadMeshEntities("../../neon/asset/dasan613_learning1/dasan613_tsdf3.obj", shader,
                                                     "../../neon/asset/dasan613_learning1/", glm::vec3(1.f), false);
    // sceneRootEntity->setEulerAngles(glm::vec3(math::toRadians(90.f), math::toRadians(0.f), math::toRadians(0.f)));
    std::cout << "min: " << sceneRootEntity->getBBox().min() << std::endl;
    std::cout << "max: " << sceneRootEntity->getBBox().max() << std::endl;

    // Point Clout Entity
    std::string pcEntityName = "PointCloud";
    std::string cloudFilename = "../../neon/asset/dasan613_learning1/cloud_learning1_4M_normals.ply";

    // Load point cloud before loadMeshEntities to give a name
    // This is ok because loadMeshEntities() first tries to find
    // whether the model given with the filename is already loaded
    ResourceManager::getModel(cloudFilename)->name = pcEntityName;
    auto pcEntityTransform =
        ECSUtil::loadMeshEntities(cloudFilename, shader, "../../neon/asset/dasan613_learning1/", glm::vec3(1.f), false);
    pcEntityTransform->getComponent<MeshRenderer>()->getMesh()->setRenderMode(GL_POINTS, 0);
    // pcEntityTransform->setEulerAngles(glm::vec3(math::toRadians(90.f), math::toRadians(0.f), math::toRadians(0.f)));
    pcEntityTransform->setPosition(glm::vec3(m_scenePosition));

    auto pcEntity = ECS::getEntityByName(pcEntityName);
    pcEntity.setActive(false);

    // Virtual sphere
    // m_sphere = EntityCreator::createSphere("virtualObject", glm::vec3(0), glm::vec3(1.f));
    std::string virtualObjectDir = "../../neon/asset/mesh/";
    std::string virtualObjectFilename = "sphere.obj";
    std::string virtualEntityName = "virtualObject";
    ResourceManager::getModel(virtualObjectDir + virtualObjectFilename)->name = virtualEntityName;
    virtualTransform = ECSUtil::loadMeshEntities(virtualObjectDir + virtualObjectFilename, shader, virtualObjectDir,
                                                 glm::vec3(5.f), false);
    virtualTransform->setPosition(glm::vec3(5.4f, 6.35f, -4.57f));

    auto sphereMaterial = EntityCreator::createMaterial();
    sphereMaterial->setFloat("u_shininess", 255.0f);
    sphereMaterial->setColor("u_color", glm::vec4(1.0f, 1.0f, 1.0f, 1.0f));
    sphereMaterial->setColor("u_emissionColor", glm::vec3(0.0f));
    sphereMaterial->setColor("u_specularColor", glm::vec3(1.f));

    auto childTransform = virtualTransform->getChildren()[0];
    childTransform->getComponent<MeshRenderer>()->setMaterial(sphereMaterial, 0);

    // Virtual should be set at the root entity
    m_sphere = ECS::getEntityByName(virtualEntityName);
    m_sphere.setVirtual(true);
    m_sphere.setActive(true);

    // // Virtual Buddha
    // auto vo1 = ResourceManager::getModel("meshes/sphere.ply");
    // vo1->name = "vo1";
    // auto vo2 = ResourceManager::getModel("meshes/buddha2.ply");
    // vo2->name = "vo2";

    // auto virtualObject = ECS::createEntity("virtualObject");
    // virtualObject.addComponent<Transform>();
    // auto parentTransform = virtualObject.getComponent<Transform>();

    // // buddha 1
    // auto voTransform1 = ECSUtil::loadMeshEntities(vo1.get(), shader, "", glm::vec3(7.5f), true);
    // // auto voTransform1 = EntityCreator::createSphere("virtualObject", glm::vec3(0),
    // // glm::vec3(1.f)).getComponent<Transform>();
    // auto buddhaMaterial1 = EntityCreator::createMaterial();
    // buddhaMaterial1->setFloat("u_shininess", 255.0);
    // buddhaMaterial1->setColor("u_color", glm::vec4(1));
    // buddhaMaterial1->setColor("u_emissionColor", glm::vec3(0.0f));
    // buddhaMaterial1->setColor("u_specularColor", glm::vec3(1.f));
    // voTransform1->getOwner().getComponent<MeshRenderer>()->setMaterial(buddhaMaterial1, 0);
    // voTransform1->getOwner().setVirtual(true);
    // voTransform1->getOwner().setActive(true);
    // voTransform1->setParent(parentTransform);

    // buddha 2
    // auto voTransform2 = ECSUtil::loadMeshEntities(vo2.get(), shader, "", glm::vec3(7.5f), true);
    // auto buddhaMaterial2 = EntityCreator::createMaterial();
    // buddhaMaterial2->setFloat("u_shininess", 255.0);
    // buddhaMaterial2->setColor("u_color", glm::vec4(1, 0, 0, 1));
    // buddhaMaterial2->setColor("u_emissionColor", glm::vec3(0.0f));
    // buddhaMaterial2->setColor("u_specularColor", glm::vec3(1.f));
    // voTransform2->getOwner().getComponent<MeshRenderer>()->setMaterial(buddhaMaterial2, 0);
    // voTransform2->getOwner().setVirtual(true);
    // voTransform2->getOwner().setActive(true);
    // voTransform2->setLocalPosition(glm::vec3(1, -0.5, 0));
    // voTransform2->setLocalEulerAngles(glm::radians(glm::vec3(180.f, 0.f, 0.f)));
    // voTransform2->setParent(parentTransform);

    // parentTransform->setPosition(glm::vec3(0.85f, 0.5f, -0.5f));
    // parentTransform->getOwner().setVirtual(true);

    // virtualTransform = parentTransform;

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
    static size_t frame = 0;
    if (DEMO_SETTINGS.animateCamera) {
        auto cameraTransform = MainCamera.getOwner().getComponent<Transform>();
        static glm::vec3 initialPos = cameraTransform->getPosition();
        static glm::quat initialRot = cameraTransform->getRotation();

        static std::shared_ptr<TransformCommand> tCommand = std::make_shared<TransformCommand>(
            cameraTransform, translations[frame], translations[frame + 1], rotations[frame], rotations[frame + 1],
            0.0333f / DEMO_SETTINGS.cameraSpeed);
        static CommandChain commandChain({tCommand}, false);
        if (commandChain.done()) {
            if (frame < rotations.size() - 1) {
                m_engine->requestScreenshot();
                m_guiEnabled = false;
                frame++;
                tCommand = std::make_shared<TransformCommand>(
                    cameraTransform, translations[frame], translations[frame + 1], rotations[frame],
                    rotations[frame + 1], 0.0333f / DEMO_SETTINGS.cameraSpeed);
                commandChain = CommandChain({tCommand}, false);
            } else {
                std::cout << "Go to initial" << std::endl;
                // Set to start position and location
                cameraTransform->setPosition(initialPos);
                cameraTransform->setRotation(initialRot);
                DEMO_SETTINGS.animateCamera.value = false;
            }
        }

        commandChain(Time::deltaTime());
    } else {
        frame = 0;
    }
}

void VoxelConeTracingDemo::updateCameraClipRegions() {
    m_clipRegionBBoxes.clear();
    auto sceneEntity = ECS::getEntityByName("name_scene");
    // glm::vec3 center = MainCamera->getPosition();
    // glm::vec3 center = sceneEntity.getComponent<Transform>()->getPosition();

    // Align with virtual object's clip region
    glm::vec3 center = virtualTransform->getPosition();

    for (size_t i = 0; i < CLIP_REGION_COUNT; ++i)
        m_clipRegionBBoxes.push_back(getBBox(i, center, m_clipRegionBBoxExtentL0));
}

void VoxelConeTracingDemo::updateVirtualClipRegions() {
    // add bbox for each clip level i around the virtual object
    m_virtualClipRegionBBoxes.clear();
    for (size_t i = 0; i < VIRTUAL_CLIP_REGION_COUNT; ++i) {
        m_virtualClipRegionBBoxes.push_back(
            getBBox(i, virtualTransform->getPosition(), m_virtualClipRegionBBoxExtentL0 / exp2f(i)));
    }
}

ClipmapUpdatePolicy::Type VoxelConeTracingDemo::getSelectedClipmapUpdatePolicyType() const {
    if (GI_SETTINGS.updateOneClipLevelPerFrame) return ClipmapUpdatePolicy::Type::ONE_PER_FRAME_PRIORITY;

    return ClipmapUpdatePolicy::Type::ALL_PER_FRAME;
}
