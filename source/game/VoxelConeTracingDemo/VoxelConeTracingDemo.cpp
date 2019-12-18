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
#include "engine/rendering/renderPasses/ForwardScenePass.h"
#include "engine/rendering/voxelConeTracing/settings/VoxelConeTracingSettings.h"
#include "engine/util/commands/RotationCommand.h"
#include "engine/util/commands/MaterialCommand.h"
#include "engine/util/commands/CommandChain.h"
#include "engine/util/QueryManager.h"
#include "engine/rendering/renderer/MeshRenderers.h"
#include "engine/util/ECSUtil/EntityCreator.h"

VoxelConeTracingDemo::VoxelConeTracingDemo()
{
    Input::subscribe(this);

    Random::randomize();
}

void VoxelConeTracingDemo::initUpdate()
{
    ResourceManager::setShaderIncludePath("shaders");        

    createDemoScene();

    DebugRenderer::init();

    init3DVoxelTextures();

    m_renderPipeline = std::make_unique<RenderPipeline>(MainCamera);
    m_gui = std::make_unique<VoxelConeTracingGUI>(m_renderPipeline.get());
    m_clipmapUpdatePolicy = std::make_unique<ClipmapUpdatePolicy>(ClipmapUpdatePolicy::Type::ONE_PER_FRAME_PRIORITY, CLIP_REGION_COUNT);

    // Set render pipeline input
    m_renderPipeline->putPtr("VoxelOpacity", &m_voxelOpacity);
    m_renderPipeline->putPtr("VoxelRadiance", &m_voxelRadiance);
    m_renderPipeline->putPtr("ClipRegionBBoxes", &m_clipRegionBBoxes);
    m_renderPipeline->putPtr("VirtualVoxelOpacity", &m_virtualVoxelOpacity);
    m_renderPipeline->putPtr("VirtualVoxelRadiance", &m_virtualVoxelRadiance);
    m_renderPipeline->putPtr("VirtualClipRegionBBoxes", &m_virtualClipRegionBBoxes);
    m_renderPipeline->putPtr("ClipmapUpdatePolicy", m_clipmapUpdatePolicy.get());
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
        std::make_shared<ForwardScenePass>());

    // RenderPass initializations
    m_renderPipeline->getRenderPass<VoxelizationPass>()->init(m_clipRegionBBoxExtentL0);

    m_initializing = false;
}

void VoxelConeTracingDemo::update()
{
#ifdef CGLAB
    static bool once = true;
#endif

    m_clipmapUpdatePolicy->setType(getSelectedClipmapUpdatePolicyType());
    m_clipmapUpdatePolicy->update();

    moveCamera(Time::deltaTime());

    glFrontFace(GL_CW);

    updateCameraClipRegions();

    for (auto c : ECS::getEntitiesWithComponents<CameraComponent>())
    {
        c.getComponent<CameraComponent>()->updateViewMatrix();
    }

    bool giPipeline = RENDERING_SETTINGS.pipeline.asString() == "GI";
    bool forwardPipeline = RENDERING_SETTINGS.pipeline.asString() == "Forward";

    if (giPipeline)
    {
        m_renderPipeline->getRenderPass<ForwardScenePass>()->setEnabled(false);
        m_renderPipeline->getRenderPass<GIPass>()->setEnabled(true);
        m_renderPipeline->getRenderPass<VoxelizationPass>()->setEnabled(true);  
#ifdef CGLAB
        // Disable theses passes because several initial updates are enough for point cloud
        m_renderPipeline->getRenderPass<RadianceInjectionPass>()->setEnabled(once);
        once = false;
#else
        m_renderPipeline->getRenderPass<VoxelizationPass>()->setEnabled(true);
        m_renderPipeline->getRenderPass<RadianceInjectionPass>()->setEnabled(true);
#endif
        m_renderPipeline->getRenderPass<SceneGeometryPass>()->setEnabled(true);
    }

    if (forwardPipeline)
    {
        m_renderPipeline->getRenderPass<ForwardScenePass>()->setEnabled(true);
        m_renderPipeline->getRenderPass<GIPass>()->setEnabled(false);
        m_renderPipeline->getRenderPass<VoxelizationPass>()->setEnabled(true);
        m_renderPipeline->getRenderPass<RadianceInjectionPass>()->setEnabled(false);
        m_renderPipeline->getRenderPass<SceneGeometryPass>()->setEnabled(false);
    }

    animateDirLight();
    animateSphereRoughness();

    GL::setViewport(MainCamera->getViewport());

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    DebugRenderInfo info(MainCamera->view(), MainCamera->proj(), MainCamera->getPosition());
    DebugRenderer::begin(info);
    DebugRenderer::update();

    glEnable(GL_DEPTH_TEST);

    m_renderPipeline->update();

    ECS::lateUpdate();

    if (m_guiEnabled)
        m_gui->update();

    m_gui->onVoxelVisualization();

    DebugRenderer::end();
}

void VoxelConeTracingDemo::moveCamera(Seconds deltaTime) const
{
    float speed = DEMO_SETTINGS.cameraSpeed * deltaTime;

    if (Input::isKeyDown(SDL_SCANCODE_W))
        MainCamera->walk(speed);

    if (Input::isKeyDown(SDL_SCANCODE_S))
        MainCamera->walk(-speed);

    if (Input::isKeyDown(SDL_SCANCODE_A))
        MainCamera->strafe(-speed);

    if (Input::isKeyDown(SDL_SCANCODE_D))
        MainCamera->strafe(speed);

    if (Input::isKeyDown(SDL_SCANCODE_Q))
        MainCamera->elevate(-speed);

    if (Input::isKeyDown(SDL_SCANCODE_E))
        MainCamera->elevate(speed);

    if (Input::rightDrag().isDragging())
    {
        float dx = -Input::rightDrag().getDragDelta().x;
        float dy = Input::rightDrag().getDragDelta().y;

        MainCamera->pitch(math::toRadians(dy * 0.1f));
        MainCamera->rotateY(math::toRadians(dx * 0.1f));
        SDL_SetRelativeMouseMode(SDL_TRUE);
    }
    else
    {
        SDL_SetRelativeMouseMode(SDL_FALSE);
    }
}

void VoxelConeTracingDemo::onKeyDown(SDL_Keycode keyCode)
{
    switch (keyCode)
    {
    case SDLK_F1:
        m_guiEnabled = !m_guiEnabled;
        break;
    case SDLK_f:
        RENDERING_SETTINGS.pipeline.curItem = 1;
        break;
    case SDLK_g:
        RENDERING_SETTINGS.pipeline.curItem = 0;
        break;
    case SDLK_F5:
        m_engine->requestScreenshot();
        break;
    default: break;
    }
}

void VoxelConeTracingDemo::init3DVoxelTextures()
{
    GLint filter = GL_LINEAR;
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

    m_voxelOpacity.create(resolutionWithBorder * FACE_COUNT, CLIP_REGION_COUNT * resolutionWithBorder, resolutionWithBorder, 
        GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE, Texture3DSettings::Custom);
    m_voxelOpacity.bind();
    m_voxelOpacity.setParameteri(GL_TEXTURE_WRAP_S, wrapS);
    m_voxelOpacity.setParameteri(GL_TEXTURE_WRAP_T, wrapT);
    m_voxelOpacity.setParameteri(GL_TEXTURE_WRAP_R, wrapR);
    m_voxelOpacity.setParameteri(GL_TEXTURE_MIN_FILTER, filter);
    m_voxelOpacity.setParameteri(GL_TEXTURE_MAG_FILTER, filter);

    m_voxelRadiance.create(resolutionWithBorder * FACE_COUNT, CLIP_REGION_COUNT * resolutionWithBorder, resolutionWithBorder, 
        GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE, Texture3DSettings::Custom);
    m_voxelRadiance.bind();
    m_voxelRadiance.setParameteri(GL_TEXTURE_WRAP_S, wrapS);
    m_voxelRadiance.setParameteri(GL_TEXTURE_WRAP_T, wrapT);
    m_voxelRadiance.setParameteri(GL_TEXTURE_WRAP_R, wrapR);
    m_voxelRadiance.setParameteri(GL_TEXTURE_MIN_FILTER, filter);
    m_voxelRadiance.setParameteri(GL_TEXTURE_MAG_FILTER, filter);

    m_virtualVoxelOpacity.create(resolutionWithBorder * FACE_COUNT, CLIP_REGION_COUNT * resolutionWithBorder, resolutionWithBorder, 
        GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE, Texture3DSettings::Custom);
    m_virtualVoxelOpacity.bind();
    m_virtualVoxelOpacity.setParameteri(GL_TEXTURE_WRAP_S, wrapS);
    m_virtualVoxelOpacity.setParameteri(GL_TEXTURE_WRAP_T, wrapT);
    m_virtualVoxelOpacity.setParameteri(GL_TEXTURE_WRAP_R, wrapR);
    m_virtualVoxelOpacity.setParameteri(GL_TEXTURE_MIN_FILTER, filter);
    m_virtualVoxelOpacity.setParameteri(GL_TEXTURE_MAG_FILTER, filter);

    m_virtualVoxelRadiance.create(resolutionWithBorder * FACE_COUNT, CLIP_REGION_COUNT * resolutionWithBorder, resolutionWithBorder, 
        GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE, Texture3DSettings::Custom);
    m_virtualVoxelRadiance.bind();
    m_virtualVoxelRadiance.setParameteri(GL_TEXTURE_WRAP_S, wrapS);
    m_virtualVoxelRadiance.setParameteri(GL_TEXTURE_WRAP_T, wrapT);
    m_virtualVoxelRadiance.setParameteri(GL_TEXTURE_WRAP_R, wrapR);
    m_virtualVoxelRadiance.setParameteri(GL_TEXTURE_MIN_FILTER, filter);
    m_virtualVoxelRadiance.setParameteri(GL_TEXTURE_MAG_FILTER, filter);

    GL_ERROR_CHECK();
}

/**
 * return: bounding box of clipmapLevel in world coordinate around center
 */
BBox VoxelConeTracingDemo::getBBox(size_t clipmapLevel, glm::vec3 center) const
{
    float halfSize = 0.5f * m_clipRegionBBoxExtentL0 * std::exp2f(float(clipmapLevel));

    return BBox(center - halfSize, center + halfSize);
}

void VoxelConeTracingDemo::createDemoScene()
{
    m_scenePosition = glm::vec3(0.0f);

    Entity camera = ECS::createEntity("Camera");
    camera.addComponent<Transform>();
    camera.addComponent<CameraComponent>();

    auto camComponent = camera.getComponent<CameraComponent>();
    auto camTransform = camera.getComponent<Transform>();

    MainCamera = camComponent;

    camComponent->setPerspective(45.0f, float(Screen::getWidth()), float(Screen::getHeight()), 0.3f, 30.0f);

#ifdef CGLAB
    // glm::vec3 cameraPositionOffset(1.f, 1.406f, -0.639f);
    glm::vec3 cameraPositionOffset(0.428, 2.249, 0.381);
    camTransform->setPosition(m_scenePosition + cameraPositionOffset);
    // camTransform->setEulerAngles(glm::vec3(math::toRadians(50.f), math::toRadians(150.f), math::toRadians(0.f)));
    camTransform->setEulerAngles(glm::vec3(math::toRadians(32.5), math::toRadians(150.f), math::toRadians(0.f)));
#else
    // For sponza
    // glm::vec3 cameraPositionOffset(8.625f, 6.593f, -0.456f);
    // camTransform->setEulerAngles(glm::vec3(math::toRadians(10.236f), math::toRadians(-66.0f), 0.0f));

    // glm::vec3 cameraPositionOffset(6.155, 5.614, 0.322);
    // camTransform->setEulerAngles(glm::vec3(math::toRadians(29.937), math::toRadians(86.204), 0.0f));

    // Ringing artifacts 
    glm::vec3 cameraPositionOffset(5.807, 4.491, -1.488);
    camTransform->setEulerAngles(glm::vec3(math::toRadians(24.738), math::toRadians(44.305), 0.0f));
    

    camTransform->setPosition(m_scenePosition + cameraPositionOffset);

#endif

    m_engine->registerCamera(camComponent);

    auto shader = ResourceManager::getShader("shaders/forwardShadingPass.vert", "shaders/forwardShadingPass.frag", {"in_pos", "in_normal", "in_tangent", "in_bitangent", "in_uv"});

    // CGLAB
#ifdef CGLAB
    ResourceManager::getModel("cglab/dasan613.obj")->name = "dasan613.obj";
    auto sceneRootEntity = ECSUtil::loadMeshEntities("cglab/dasan613.obj", shader, "cglab/", glm::vec3(1.f), true);
    sceneRootEntity->setEulerAngles(glm::vec3(math::toRadians(90.f), math::toRadians(0.f), math::toRadians(0.f)));

    // Point Clout Entity
    std::string pcEntityName = "PointCloud";
    std::string cloudFilename = "cglab/cloud_binary.ply";

    // Load point cloud before loadMeshEntities to give a name
    // This is ok because loadMeshEntities() first tries to find 
    // whether the model given with the filename is already loaded
    ResourceManager::getModel(cloudFilename)->name = pcEntityName; 
    auto pcEntityTransform = ECSUtil::loadMeshEntities(cloudFilename, shader, "cglab/", glm::vec3(1.f), false);
    
    pcEntityTransform->getComponent<MeshRenderer>()->getMesh()->setRenderMode(GL_POINTS, 0);
    pcEntityTransform->setEulerAngles(glm::vec3(math::toRadians(90.f), math::toRadians(0.f), math::toRadians(0.f)));
    pcEntityTransform->setPosition(glm::vec3(m_scenePosition));

    auto pcEntity = ECS::getEntityByName(pcEntityName);
    pcEntity.setActive(false);
#else
    // For sponza
    auto sceneRootEntity = ECSUtil::loadMeshEntities("meshes/sponza_obj/sponza.obj", shader, "textures/sponza_textures/", glm::vec3(0.01f), true);
#endif

    // // Virtual sphere
    // m_sphere = EntityCreator::createSphere("sphere", glm::vec3(0), glm::vec3(1.f));
    // auto sphereMaterial = EntityCreator::createMaterial();
    // sphereMaterial->setFloat("u_shininess", 255.0f);
    // sphereMaterial->setColor("u_color", glm::vec4(0.0f, 0.0f, 0.0f, 1.0f));
    // sphereMaterial->setColor("u_emissionColor", glm::vec3(0.0f));
    // sphereMaterial->setColor("u_specularColor", glm::vec3(1.f));
    // m_sphere.getComponent<MeshRenderer>()->setMaterial(sphereMaterial, 0);

    // Virtual Buddha
    ResourceManager::getModel("meshes/buddha/buddha.ply")->name = "buddha"; 
    buddhaTransform = ECSUtil::loadMeshEntities("meshes/buddha/buddha.ply", shader, "", glm::vec3(10.f), true);
    auto buddhaMaterial = EntityCreator::createMaterial();
    buddhaMaterial->setFloat("u_shininess", 255.0f);
    buddhaMaterial->setColor("u_color", glm::vec4(0.0f, 0.0f, 0.0f, 1.0f));
    buddhaMaterial->setColor("u_emissionColor", glm::vec3(0.0f));
    buddhaMaterial->setColor("u_specularColor", glm::vec3(1.0f));
    buddhaTransform->getOwner().getComponent<MeshRenderer>()->setMaterial(buddhaMaterial, 0);
    buddhaTransform->getOwner().setVirtual(true);

#ifdef CGLAB
    // m_sphere.getComponent<Transform>()->setPosition(glm::vec3(1.35, 0.45, -1.3));
    // m_sphere.getComponent<Transform>()->setLocalScale(glm::vec3(0.7f));
    buddhaTransform->setPosition(glm::vec3(1.35, 0.95, -1.3));
    buddhaTransform->setEulerAngles(glm::vec3(0.f, 90.f, 0.f));
#else
    // Original
    // m_sphere.getComponent<Transform>()->setPosition(glm::vec3(7.035, 5.092, 0.396));

    // Ringing artifacts
    // m_sphere.getComponent<Transform>()->setPosition(glm::vec3(6.485, 3.942, -0.554));

    // Buddha
    buddhaTransform->setPosition(glm::vec3(7.035, 5.092, 0.396));
    buddhaTransform->setEulerAngles(glm::vec3(0.f, 90.f, 0.f));
#endif

    if (sceneRootEntity)
        sceneRootEntity->setPosition(glm::vec3(m_scenePosition));

#ifndef CGLAB
    m_directionalLight = ECS::createEntity("Directional Light");
    m_directionalLight.addComponent<DirectionalLight>();
    m_directionalLight.addComponent<Transform>();
    m_directionalLight.getComponent<Transform>()->setPosition(glm::vec3(0.0f, 20.0f, 0.f) + m_scenePosition);
    m_directionalLight.getComponent<Transform>()->setEulerAngles(glm::vec3(math::toRadians(111.5f), math::toRadians(3.5f), 0.f));
    m_directionalLight.getComponent<DirectionalLight>()->intensity = 1.5f;
#endif
}

void VoxelConeTracingDemo::animateDirLight()
{
    if (m_directionalLight && DEMO_SETTINGS.animateLight)
    {
        auto lightTransform = m_directionalLight.getComponent<Transform>();
        glm::vec3 lightPos = lightTransform->getPosition();
        glm::quat startRotation = glm::quat(glm::transpose(glm::lookAt(lightPos, glm::vec3(0.0f, 0.0f, -10.0f), glm::vec3(0.0f, 1.0f, 0.0f))));
        glm::quat dstRotation = glm::quat(glm::transpose(glm::lookAt(lightPos, glm::vec3(0.0f, 0.0f, 10.0f), glm::vec3(0.0f, 1.0f, 0.0f))));

        static std::shared_ptr<RotationCommand> rotationCommand = std::make_shared<RotationCommand>(lightTransform, startRotation, dstRotation, 10.0f);
        static std::shared_ptr<RotationCommand> rotationCommand2 = std::make_shared<RotationCommand>(lightTransform, dstRotation, startRotation, 10.0f);
        static CommandChain commandChain({rotationCommand, rotationCommand2}, true);

        commandChain(Time::deltaTime());
    }
}

void VoxelConeTracingDemo::animateSphereRoughness()
{
    if (m_sphere && DEMO_SETTINGS.animateSphere)
    {
        auto sphereRenderer = m_sphere.getComponent<MeshRenderer>();
        auto material = sphereRenderer->getMaterial(0);
        auto startRoughness = glm::vec4(255.f);
        auto dstRoughness = glm::vec4(0);

        static std::shared_ptr<MaterialCommand> materialCommand = std::make_shared<MaterialCommand>(sphereRenderer, "u_shininess", startRoughness, dstRoughness, 5.0f);
        static std::shared_ptr<MaterialCommand> materialCommand2 = std::make_shared<MaterialCommand>(sphereRenderer, "u_shininess", dstRoughness, startRoughness, 5.0f);
        static CommandChain commandChain({materialCommand, materialCommand2}, true);

        commandChain(Time::deltaTime());
    }
}

void VoxelConeTracingDemo::updateCameraClipRegions()
{
    m_clipRegionBBoxes.clear();
    glm::vec3 center = MainCamera->getPosition();
    for (size_t i = 0; i < CLIP_REGION_COUNT; ++i)
        m_clipRegionBBoxes.push_back(getBBox(i, center));
}

void VoxelConeTracingDemo::updateVirtualClipRegions()
{
    // add bbox for each clip level i around the virtual object
    m_virtualClipRegionBBoxes.clear();
    for (size_t i = 0; i < CLIP_REGION_COUNT; ++i) {
        // m_virtualClipRegionBBoxes.push_back(buddhaTransform->getBBox());
        m_virtualClipRegionBBoxes.push_back(getBBox(i, buddhaTransform->getPosition()));
    }
}

ClipmapUpdatePolicy::Type VoxelConeTracingDemo::getSelectedClipmapUpdatePolicyType() const
{
    if (GI_SETTINGS.updateOneClipLevelPerFrame)
        return ClipmapUpdatePolicy::Type::ONE_PER_FRAME_PRIORITY;

    return ClipmapUpdatePolicy::Type::ALL_PER_FRAME;
}
