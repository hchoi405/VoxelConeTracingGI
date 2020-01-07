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

const float MAX_TRACE_DISTANCE = 30.0;
const float MIN_STEP_FACTOR = 0.2;
const float MIN_SPECULAR_APERTURE = 0.05;

uniform sampler3D u_voxelRadiance;
uniform sampler3D u_voxelOpacity;
uniform sampler2D u_depthTexture;
uniform sampler2D u_normalMap;

uniform uint u_volumeDimension;
uniform float u_voxelSizeL0;
uniform vec3 u_volumeCenterL0;
uniform float u_traceStartOffset;
uniform float u_occlusionDecay = 1.0;
uniform float u_stepFactor;
uniform float u_viewAperture;
uniform float u_hitpointOffset;

uniform vec3 u_eyePos;
uniform mat4 u_viewProjInv;

layout (location = 0) out vec4 out_color;

float getMinLevel(vec3 posW)
{
    float distanceToCenter = length(u_volumeCenterL0 - posW);
    float minRadius = u_voxelSizeL0 * u_volumeDimension * 0.5;
    float minLevel = log2(distanceToCenter / minRadius);  
    minLevel = max(0.0, minLevel);
    
    float radius = minRadius * exp2(ceil(minLevel));
    float f = distanceToCenter / radius;
    
    // Smoothly transition from current level to the next level
    float transitionStart = 0.5;
    float c = 1.0 / (1.0 - transitionStart);
    
    return f > transitionStart ? ceil(minLevel) + (f - transitionStart) * c : ceil(minLevel);
}

vec3 worldToSamplePos(vec3 posW, float extent)
{
	return (fract(posW / extent) * u_volumeDimension + vec3(BORDER_WIDTH)) / (float(u_volumeDimension) + 2.0 * BORDER_WIDTH);
}

float dx[8] = {0, 0, 0, 0,  1,  1,  1,  1};
float dy[8] = {0, 0,  1,  1, 0, 0,  1,  1};
float dz[8] = {0,  1, 0,  1, 0,  1, 0,  1};

// All input texture should have GL_NEAREST filter
vec4 sampleClipmapTexture2(sampler3D radianceTexture, sampler3D opacityTexture, vec3 posW, int clipmapLevel, vec3 faceOffsets, vec3 weight)
{
    float voxelSize = u_voxelSizeL0 * exp2(clipmapLevel);
    float extent = voxelSize * u_volumeDimension;

    vec3 neighborPos[8];
    vec3 neighborColor[8];
    bool neighborExist[8];

    // Trilinearly interpolate
    vec3 sum = vec3(0.f);

    vec3 a = round(posW / voxelSize) * voxelSize;   // nearest voxel
    for (uint i=0; i<8; ++i)
    {
        neighborExist[i] = false;
        neighborPos[i] = posW + vec3(dx[i], dy[i], dz[i]) * voxelSize;
        neighborColor[i] = vec3(0.0);

#ifdef VOXEL_TEXTURE_WITH_BORDER
        neighborPos[i] = (fract(neighborPos[i] / extent) * u_volumeDimension + vec3(BORDER_WIDTH)) / (float(u_volumeDimension) + 2.0 * BORDER_WIDTH);
#else
        neighborPos[i] = fract(neighborPos[i] / extent);
#endif 
        neighborPos[i].y += clipmapLevel;
        neighborPos[i].y *= CLIP_LEVEL_COUNT_INV;
        neighborPos[i].x *= FACE_COUNT_INV;

        float opaque[3] = {
            texture(opacityTexture, neighborPos[i] + vec3(faceOffsets.x, 0.0, 0.0)).r,
            texture(opacityTexture, neighborPos[i] + vec3(faceOffsets.y, 0.0, 0.0)).r,
            texture(opacityTexture, neighborPos[i] + vec3(faceOffsets.z, 0.0, 0.0)).r
        };

        uint faceCount = 0;
        for (uint j=0; j<3; ++j)
        {
            if (opaque[j] > 0.0) {
                neighborColor[i] += texture(radianceTexture, neighborPos[i] + vec3(faceOffsets[j], 0.0, 0.0)).rgb;
                faceCount++;
            }
        }
        if (faceCount > 0) {
            neighborColor[i] /= faceCount;
            neighborExist[i] = true;
        }
    }

    vec3 frac = fract(posW / voxelSize);
    
    bool bv01, bv23, bv0123, bv45, bv67, bv4567, bv01234567;
    vec3 v01, v23, v0123, v45, v67, v4567, v01234567;

    bv01 = true;
    if (neighborExist[0] && neighborExist[1]) {
        v01 = mix(neighborColor[0], neighborColor[1], frac.z);
    } else if (neighborExist[0]) {
        v01 = neighborColor[0];
    } else if (neighborExist[1]) {
        v01 = neighborColor[1];
    } else {
        bv01 = false;
    }

    bv23 = true;
    if (neighborExist[2] && neighborExist[3]) {
        v23 = mix(neighborColor[2], neighborColor[3], frac.z);
    } else if (neighborExist[2]) {
        v23 = neighborColor[2];
    } else if (neighborExist[3]) {
        v23 = neighborColor[3];
    } else {
        bv23 = false;
    }

    bv0123 = true;
    if (bv01 && bv23) {
        v0123 = mix(v01, v23, frac.y);
    } else if (bv01) {
        v0123 = v01;
    } else if (bv23) {
        v0123 = v23;
    } else {
        bv0123 = false;
    }

    bv45 = true;
    if (neighborExist[4] && neighborExist[5]) {
        v45 = mix(neighborColor[4], neighborColor[5], frac.z);
    } else if (neighborExist[4]) {
        v45 = neighborColor[4];
    } else if (neighborExist[5]) {
        v45 = neighborColor[5];
    } else {
        bv45 = false;
    }

    bv67 = true;
    if (neighborExist[6] && neighborExist[7]) {
        v67 = mix(neighborColor[6], neighborColor[7], frac.z);
    } else if (neighborExist[6]) {
        v67 = neighborColor[6];
    } else if (neighborExist[7]) {
        v67 = neighborColor[7];
    } else {
        bv67 = false;
    }

    bv4567 = true;
    if (bv45 && bv67) {
        v4567 = mix(v45, v67, frac.y);
    } else if (bv45) {
        v4567 = v45;
    } else if (bv67) {
        v4567 = v67;
    } else {
        bv4567 = false;
    }

    if (bv0123 && bv4567) {
        v01234567 = mix(v0123, v4567, frac.x);
    } else if (bv0123) {
        v01234567 = v0123;
    } else if (bv4567) {
        v01234567 = v4567;
    } else {
        v01234567 = vec3(0.0);
    }

    sum = v01234567;
    
    return clamp(vec4(sum, 1.0), 0.0, 1.0);
}

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
    
    vec4 lowSample = sampleClipmapTexture2(clipmapTexture, u_voxelOpacity, posW, lowerLevel, faceOffsets, weight);
    
	if (lowerLevel == upperLevel)
        return lowSample;
	
    vec4 highSample = sampleClipmapTexture2(clipmapTexture, u_voxelOpacity, posW, upperLevel, faceOffsets, weight);
	
    return mix(lowSample, highSample, fract(curLevel));
}

vec3 worldPosFromDepth(float depth)
{
    vec4 p = vec4(In.texCoords, depth, 1.0);
    p.xyz = p.xyz * 2.0 - 1.0;
    p = u_viewProjInv * p;
    return p.xyz / p.w;
}

void main() 
{
    float depth = texture2D(u_depthTexture, In.texCoords).r;
    if (depth == 1.0)
        discard;
    vec3 posW = worldPosFromDepth(depth);
    vec3 normal = texture(u_normalMap, In.texCoords).rgb;
    posW += u_hitpointOffset * normalize(normal) * u_voxelSizeL0;
    vec3 view = normalize(posW - u_eyePos);
    // vec3 view = normalize(normal);
    ivec3 faceIndices = computeVoxelFaceIndices(view);

    // Real weight
    vec3 w = view;

    // We can use this becase all faces has the same value (because of storeVoxelColorAtomicRGBA8Avg6Faces)
    // vec3 w = normalize(vec3(1.0));

    // Weired weight normalization?
    w = w*w;
    out_color = sampleClipmapLinearly(u_voxelRadiance, posW, 0, faceIndices, w);
    // out_color = vec4(castCone(u_eyePos, view, u_viewAperture, MAX_TRACE_DISTANCE, 0).rgb, 1.0);
    // indirectContribution += castCone(startPos, DIFFUSE_CONE_DIRECTIONS[i], DIFFUSE_CONE_APERTURE ,MAX_TRACE_DISTANCE, minLevel) * cosTheta;
}
