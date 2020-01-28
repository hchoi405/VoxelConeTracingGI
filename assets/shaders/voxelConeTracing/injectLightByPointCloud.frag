#version 430
#extension GL_ARB_shading_language_include : enable
#extension GL_ARB_shader_image_load_store : require
#include "/voxelConeTracing/voxelizationFrag.glsl"
#include "/voxelConeTracing/common.glsl"
#include "/shadows/shadows.glsl"
#include "/voxelConeTracing/conversion.glsl"

in Geometry
{
    vec3 posW;
    vec3 normalW;
    vec2 uv;
    vec3 color;
} In;

uniform sampler2D u_diffuseTexture0;
uniform sampler2D u_emissionMap0;
uniform sampler2D u_opacityMap0;
uniform float u_hasEmissionMap;
uniform float u_hasOpacityMap;
uniform float u_hasDiffuseTexture;

uniform float u_shininess;
uniform vec4 u_color;
uniform vec3 u_emissionColor;
uniform vec3 u_specularColor;

uniform int u_normalOnly = 0;
uniform uint u_noSpecular = 0;

uniform layout(r32ui) volatile uimage3D u_voxelRadiance;
uniform layout(r32ui) volatile uimage3D u_voxelNormal;
uniform layout(r32ui) volatile uimage3D u_voxelDiffuse;
uniform layout(r32ui) volatile uimage3D u_voxelSpecularA;

void main() 
{
    vec3 posW = In.posW;
        
    if (isOutsideVoxelizationRegion(posW) || isInsideDownsampleRegion(posW))
        discard;
        
	if (u_hasOpacityMap > 0.0 && texture(u_opacityMap0, In.uv).r < 0.1)
        discard;

    vec3 normal = normalize(In.normalW);
    ivec3 faceIndices = computeVoxelFaceIndices(-normal);

    // Color
    if (u_normalOnly != 1) {
        // storeVoxelColorAtomicRGBA8Avg6Faces(u_voxelRadiance, posW, vec4(In.color, 1.0));
        storeVoxelColorAtomicRGBA8Avg(u_voxelRadiance, posW, vec4(In.color, 1.0), faceIndices, abs(normal));
        // storeVoxelColorAtomicRGBA8Avg(u_voxelRadiance, posW, vec4(1,0,0, 1.0), faceIndices, abs(normal));
    }

    // Normal
    storeVoxelColorAtomicRGBA8Avg(u_voxelNormal, posW, vec4(packNormal(normal), 1.0), faceIndices, abs(normal));
    
    // Diffuse
    storeVoxelColorAtomicRGBA8Avg(u_voxelDiffuse, posW, u_color, faceIndices, abs(normal));

    // Specular
    if (u_noSpecular == 0)
        storeVoxelColorAtomicRGBA8Avg(u_voxelSpecularA, posW, vec4(u_specularColor, packShininess(u_shininess)), faceIndices, abs(normal));

}
