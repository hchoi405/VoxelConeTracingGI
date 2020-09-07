#pragma once
#include <glm/glm.hpp>
#include <vector>
#include "engine/gui/GUIElements.h"
#include "engine/rendering/voxelConeTracing/Globals.h"

struct VCTSettings
{
    std::vector<GUIElement*> guiElements;
};

struct ShadowSettings : VCTSettings
{
    ShadowSettings() { guiElements.insert(guiElements.end(), {&usePoissonFilter, &depthBias, &radianceVoxelizationPCFRadius}); }

    CheckBox usePoissonFilter{"Use Poisson Filter", true};

    uint32_t shadowMapResolution{4096}; // 8184

    SliderFloat depthBias{"Depth Bias", 0.013f, 0.0001f, 0.1f, "%.6f"};
    SliderFloat radianceVoxelizationPCFRadius{ "Radiance Voxelization PCF Radius", 0.5f / VOXEL_RESOLUTION, 0.0f, 2.0f / VOXEL_RESOLUTION };
};

struct RenderingSettings : VCTSettings
{
    RenderingSettings() { guiElements.insert(guiElements.end(), {&wireFrame, &pipeline, &cullBackFaces, &brdfMode }); }

    CheckBox wireFrame{"Wireframe", false};
    ComboBox pipeline = ComboBox("Pipeline", {"GI", "Forward"}, 0);
    CheckBox cullBackFaces{ "Cull Back Faces", false };
    ComboBox brdfMode = ComboBox("BRDF", { "Blinn-Phong", "Cook-Torrance"}, 1);
};

struct VisualizationSettings : VCTSettings
{
    VisualizationSettings() { guiElements.insert(guiElements.end(), {&voxelVisualizationAlpha, &borderWidth, &borderColor}); }

    SliderFloat voxelVisualizationAlpha{"Voxel Alpha", 1.0f, 0.0f, 1.0f};
    SliderFloat borderWidth{"Border Width", 0.05f, 0.0f, 1.0f};
    ColorSelection borderColor{"Border Color", glm::vec4(0.5f, 0.5f, 0.5f, 1.0f)};
};

struct GISettings : VCTSettings
{
    GISettings()
    {
        guiElements.insert(guiElements.end(), {&occlusionDecay, &ambientOcclusionFactor, &stepFactor, &virtualStepFactor,
                          &realIndirectDiffuseIntensity, &virtualIndirectDiffuseIntensity, &indirectSpecularIntensity, 
                          &traceStartOffset, &traceDirectionOffset, &directLighting, &indirectDiffuseLighting, 
                          &indirectSpecularLighting, &ambientOcclusion, &radianceInjectionMode, 
                          &visualizeMinLevelSelection, &downsampleTransitionRegionSize, &updateOneClipLevelPerFrame });
    }

    // Dasan613, mirror sphere: 0.001
    // Dasan106, cube diffuse: 5
    SliderFloat occlusionDecay{"Occlusion Decay",5.f, 0.001f, 80.0f};
    SliderFloat ambientOcclusionFactor{ "Ambient Occlusion Factor", 2.0f, 0.1f, 4.0f };
    // 1.0: Dasan613
    // Dasan613, mirror sphere: 0.785
    SliderFloat stepFactor{"Step Factor", 0.785f, 0.2f, 2.0f};
    SliderFloat virtualStepFactor{"Virtual Step Factor", 0.1f, 0.1f, 1.0f};
    SliderFloat realIndirectDiffuseIntensity{"Real Indirect Diffuse Intensity", 5.f, 0.1f, 15.0f};
    // Dasan106, diffuse cube (nodirect): 4.f
    // Dasan106, diffuse cube (direct): 6.782f
    // Dasan106, diffuse cube (direct, rotate): 4.4f
    // Rise106, diffuse buddha (direct, norotate): 1.9f
    SliderFloat virtualIndirectDiffuseIntensity{"Virtual Indirect Diffuse Intensity", 1.9f, 0.1f, 15.0f};
    // Dasan613, mirror lucy: 1.34
    SliderFloat indirectSpecularIntensity{ "Indirect Specular Intensity", 1.34f, 0.1f, 3.0f };
    SliderFloat traceStartOffset{"Trace Start Offset", 1.5f, 0.0f, 8.0f};
    SliderFloat traceDirectionOffset{"Trace Direction Offset", 1.5f, 0.0f, 8.0f};
    
    CheckBox directLighting{ "Direct Lighting", true };
    CheckBox indirectDiffuseLighting{ "Indirect Diffuse Lighting", true };
    CheckBox indirectSpecularLighting{ "Indirect Specular Lighting", true };
    CheckBox ambientOcclusion{ "Ambient Occlusion", false };
    ComboBox radianceInjectionMode = ComboBox("Radiance Injection Mode", { "Conservative", "MSAA", "Point Cloud" }, 2);
    CheckBox visualizeMinLevelSelection{"Visualize Min Level Selection", false};
    SliderInt downsampleTransitionRegionSize{ "Downsample Transition Region Size", 10, 1, VOXEL_RESOLUTION / 4 };
    CheckBox updateOneClipLevelPerFrame{ "Update One Clip Level Per Frame", false };
};

struct DebugSettings : VCTSettings
{
    DebugSettings()
    {
        guiElements.insert(
            guiElements.end(),
            {&debugFlag, /* &toggleViewCone, */ &viewAperture, /*  &hitpointOffset, &virtualSelfOcclusion,
             &indirectSpecularShadow, &indirectDiffuseShadow, &irradianceOnly, */ &secondBounce, &secondIndirectDiffuse,
             &secondIndirectSpecular, &realReflectance, &renderReal, &renderVirtual, /* &ambientSecondIntensity,
             &extraStep, */ &glassEta, &phongShininess, &rotateCone, &localRatio, &exclueEmptyFace, &materialMode,
             &subsample});
    }

    CheckBox toggleViewCone{"Toggle view-based cone tracing", false};
    SliderFloat viewAperture{"Apertuer of View Cone", 0.05f, 0.0f, 1.0f};
    SliderFloat hitpointOffset{"Offset of hitpoint to normal direction", 10.f, 0.0f, 10.0f};
    CheckBox virtualSelfOcclusion{"Self occlusion for virtual object", true};
    SliderFloat indirectSpecularShadow{"indirectSpecularShadow", 0.45f, 0.01f, 10.f};
    SliderFloat indirectDiffuseShadow{"indirectDiffuseShadow", 2.5f, 0.01f, 10.f};
    CheckBox irradianceOnly{"Show irradiance (denominator) only instead of reflectance", false};
    // Dasan106, diffuse buddha: 2..f
    // Rise103, diffuse buddha: 1.f
    SliderFloat secondIndirectDiffuse{"Second bounce diffuse intenfsity", 1.f, 0.01f, 3.f};
    SliderFloat secondIndirectSpecular{"Second bounce specular intensity", 1.f, 0.01f, 3.f};
    CheckBox secondBounce{"Trace second bounce", true};
    CheckBox realReflectance{"Apply relfectance of real object", true};
    CheckBox debugFlag{"Flag for Debug", false};
    CheckBox renderReal{"Toggle whether render real fragment", true};
    CheckBox renderVirtual{"Toggle whether render virtual fragment", true};
    SliderFloat ambientSecondIntensity{"Ambient light intensity for second bounce", 1.f, 0.01f, 1.f};
    SliderInt extraStep{"Number of steps to sample after primary", 3, 0, 100};
    SliderInt subsample{"Number of subsample from primary", 1, 1, 128};
    SliderFloat glassEta{"Eta for glass", 1.5f, 1.f, 2.f};
    SliderFloat phongShininess{"Shininess of phong", 27.8974f, 0.f, 255.f};
    CheckBox rotateCone{"Rotate cone randomly", false};
    SliderFloat localRatio{"Ratio of local real", 0.5f, 0.f, 2.f};
    CheckBox exclueEmptyFace{"Exclude empty faces when rendering (custom interpol)", false};
    ComboBox materialMode = ComboBox("Material", { "Mirror", "Glass", "Diffuse", "Phong"}, 0);
};

struct DemoSettings : VCTSettings
{
    DemoSettings()
    {
        guiElements.insert(guiElements.end(), {&animateLight, &animateSphere, &animateCamera, &cameraSpeed, &animateFrame });
    }

    CheckBox animateLight{ "Animate Light", false };
    CheckBox animateSphere{ "Animate Sphere Roughness", false };
    CheckBox animateCamera{ "Animate Camera Transform", false };
    SliderFloat cameraSpeed{ "Camera Speed", 3.0f, 1.0f, 15.0f };
    // SliderInt animateFrame{ "Frame", 0, 0, 450 };
    SliderInt animateFrame{"Frame", 120, 120, 120};
};
