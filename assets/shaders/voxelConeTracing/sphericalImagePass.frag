#version 430
#extension GL_ARB_shading_language_include : enable

#include "/voxelConeTracing/settings.glsl"
#include "/BRDF.glsl"
#include "/voxelConeTracing/common.glsl"
#include "/shadows/shadows.glsl"
#include "/intersection.glsl"
#include "/voxelConeTracing/conversion.glsl"

in Vertex
{
    vec2 texCoords;
} In;

uniform sampler2D u_positionMap;
uniform sampler3D u_voxelRadiance;

uniform uint u_volumeDimension;
uniform float u_voxelSizeL0;
uniform vec3 u_sphericalCenter;

layout (location = 0) out vec3 out_color;

vec4 sampleClipmapTexture(sampler3D clipmapTexture, vec3 posW, int clipmapLevel, vec3 faceOffsets, vec3 weight)
{
	float voxelSize = u_voxelSizeL0 * exp2(clipmapLevel);
    float extent = voxelSize * u_volumeDimension;
	
#ifdef VOXEL_TEXTURE_WITH_BORDER
	vec3 samplePos = (fract(posW / extent) * u_volumeDimension + vec3(BORDER_WIDTH)) / (float(u_volumeDimension) + 2.0 * BORDER_WIDTH);
#else
    vec3 samplePos = fract(posW / extent);
#endif

    samplePos.y += clipmapLevel;
    samplePos.y *= CLIP_LEVEL_COUNT_INV;
    samplePos.x *= FACE_COUNT_INV;

    return clamp(texture(clipmapTexture, samplePos + vec3(faceOffsets.x, 0.0, 0.0)) * weight.x +
                 texture(clipmapTexture, samplePos + vec3(faceOffsets.y, 0.0, 0.0)) * weight.y +
                 texture(clipmapTexture, samplePos + vec3(faceOffsets.z, 0.0, 0.0)) * weight.z, 0.0, 1.0);
}

vec4 sampleClipmapLinearly(sampler3D clipmapTexture, vec3 posW, float curLevel, ivec3 faceIndices, vec3 weight)
{    
    int lowerLevel = int(floor(curLevel));
    int upperLevel = int(ceil(curLevel));
    
    vec3 faceOffsets = vec3(faceIndices) * FACE_COUNT_INV;
    
    vec4 lowSample = sampleClipmapTexture(clipmapTexture, posW, lowerLevel, faceOffsets, weight);
    
	if (lowerLevel == upperLevel)
        return lowSample;
	
    vec4 highSample = sampleClipmapTexture(clipmapTexture, posW, upperLevel, faceOffsets, weight);
	
    return mix(lowSample, highSample, fract(curLevel));
}

void main() 
{
    vec3 hitPoint = texture(u_positionMap, In.texCoords).rgb;
    vec3 view = hitPoint - u_sphericalCenter;
    ivec3 faceIndices = computeVoxelFaceIndices(view);
    // vec3 w = view * view;
    vec3 w = vec3(1.0);
    out_color = sampleClipmapLinearly(u_voxelRadiance, hitPoint, 0, faceIndices, w).rgb;
}
