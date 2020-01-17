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
        guiElements.insert(guiElements.end(), {&occlusionDecay, &ambientOcclusionFactor, &stepFactor,
                          &indirectDiffuseIntensity, &indirectSpecularIntensity, &traceStartOffset,
                          &directLighting, &indirectDiffuseLighting, &indirectSpecularLighting, &ambientOcclusion,
                          &radianceInjectionMode, &visualizeMinLevelSelection, &downsampleTransitionRegionSize,
                          &updateOneClipLevelPerFrame });
    }

    SliderFloat occlusionDecay{"Occlusion Decay", 5.0f, 0.001f, 80.0f};
    SliderFloat ambientOcclusionFactor{ "Ambient Occlusion Factor", 2.0f, 0.1f, 4.0f };
    SliderFloat stepFactor{"Step Factor", 0.2f, 0.2f, 2.0f};
    SliderFloat indirectDiffuseIntensity{"Indirect Diffuse Intensity", 11.0f, 1.0f, 30.0f};
    SliderFloat indirectSpecularIntensity{ "Indirect Specular Intensity", 2.0f, 1.0f, 16.0f };
    SliderFloat traceStartOffset{"Trace Start Offset", 1.5f, 0.0f, 8.0f};
    
    CheckBox directLighting{ "Direct Lighting", true };
    CheckBox indirectDiffuseLighting{ "Indirect Diffuse Lighting", true };
    CheckBox indirectSpecularLighting{ "Indirect Specular Lighting", true };
    CheckBox ambientOcclusion{ "Ambient Occlusion", true };
    ComboBox radianceInjectionMode = ComboBox("Radiance Injection Mode", { "Conservative", "MSAA", "Point Cloud" }, 2);
    CheckBox visualizeMinLevelSelection{"Visualize Min Level Selection", false};
    SliderInt downsampleTransitionRegionSize{ "Downsample Transition Region Size", 10, 1, VOXEL_RESOLUTION / 4 };
    CheckBox updateOneClipLevelPerFrame{ "Update One Clip Level Per Frame", false };
};

struct DebugSettings : VCTSettings
{
    DebugSettings()
    {
        guiElements.insert(guiElements.end(), {
            /* &viewAperture, &hitpointOffset, &raymarchingCounter, &virtualStepFactor,
        &indirectVirtualRadius, &opacityCorrection, &counterBreak,  */
        &virtualSelfOcclusion, &indirectSpecularShadow,
        &indirectDiffuseShadow, &irradianceOnly, &secondBounce, &secondIndirectDiffuse, &realReflectance});
    }

    SliderFloat viewAperture{"Apertuer of View Cone", 0.05f, 0.0f, 1.0f};
    SliderFloat hitpointOffset{"Offset of hitpoint to normal direction", 0.f, -1.0f, 1.0f};
    SliderFloat raymarchingCounter{"Ray Marching Counter (darker part = small)", 1, 0, 10};
    SliderFloat virtualStepFactor{"virtualStepFactor", 0.2, 0.01, 10};
    SliderFloat indirectVirtualRadius{"indirectVirtualRadius", 1, 1, 256};
    CheckBox opacityCorrection{"Opacity Correction", true};
    SliderInt counterBreak{"Counter Break", 1000, 0, 1000};
    CheckBox virtualSelfOcclusion{"Self occlusion for virtual object", true};
    SliderFloat indirectSpecularShadow{"indirectSpecularShadow", 0.45f, 0.01f, 10.f};
    SliderFloat indirectDiffuseShadow{"indirectDiffuseShadow", 0.25f, 0.01f, 10.f};
    CheckBox irradianceOnly{"Show irradiance (denominator) only instead of reflectance", false};
    SliderFloat secondIndirectDiffuse{"Second bounce diffuse factor", 30.f, 1.f, 30.f};
    CheckBox secondBounce{"Trace second bounce", false};
    CheckBox realReflectance{"Apply relfectance of real object", true};
};

struct DemoSettings : VCTSettings
{
    DemoSettings()
    {
        guiElements.insert(guiElements.end(), {&animateLight, &animateSphere, /*animateCamera,*/ &cameraSpeed });
    }

    CheckBox animateLight{ "Animate Light", false };
    CheckBox animateSphere{ "Animate Sphere Roughness", false };
    CheckBox animateCamera{ "Animate Camera Transform", false };
    SliderFloat cameraSpeed{ "Camera Speed", 3.0f, 1.0f, 15.0f };
};
