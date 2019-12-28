#version 430
#extension GL_ARB_shading_language_include : enable
#extension GL_ARB_shader_image_load_store : require
#include "/voxelConeTracing/voxelizationFrag.glsl"
#include "/voxelConeTracing/common.glsl"
#include "/shadows/shadows.glsl"

in Geometry
{
    vec3 posW;
    vec3 normalW;
    vec2 uv;
    vec3 color;
} In;

uniform DirectionalLight u_directionalLights[MAX_DIR_LIGHT_COUNT];
uniform DirectionalLightShadowDesc u_directionalLightShadowDescs[MAX_DIR_LIGHT_COUNT];
uniform sampler2D u_shadowMaps[MAX_DIR_LIGHT_COUNT];
uniform int u_numActiveDirLights;
uniform float u_depthBias;
uniform float u_usePoissonFilter;

uniform sampler2D u_diffuseTexture0;
uniform sampler2D u_emissionMap0;
uniform sampler2D u_opacityMap0;
uniform float u_hasEmissionMap;
uniform float u_hasOpacityMap;
uniform float u_hasDiffuseTexture;
uniform vec4 u_color;
uniform vec3 u_emissionColor;

uniform layout(r32ui) volatile uimage3D u_voxelRadiance;

void main() 
{
    vec3 posW = In.posW;
        
    if (isOutsideVoxelizationRegion(posW) || isInsideDownsampleRegion(posW))
        discard;
        
	if (u_hasOpacityMap > 0.0 && texture(u_opacityMap0, In.uv).r < 0.1)
        discard;

    vec3 normal = normalize(In.normalW);
    ivec3 faceIndices = computeVoxelFaceIndices(-normal);
    // storeVoxelColorAtomicRGBA8Avg6Faces(u_voxelRadiance, posW, vec4(In.color, 1.0));
    storeVoxelColorAtomicRGBA8Avg(u_voxelRadiance, posW, vec4(In.color, 1.0), faceIndices, abs(normal));
}
