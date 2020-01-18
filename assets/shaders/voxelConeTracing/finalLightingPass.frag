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

struct Intersection {
    bool hit;
    vec3 position;
    vec3 normal;
    float level;
    bool hitOnly;
};

#define DIRECT_LIGHTING_BIT 1
#define INDIRECT_DIFFUSE_LIGHTING_BIT 2
#define INDIRECT_SPECULAR_LIGHTING_BIT 4
#define AMBIENT_OCCLUSION_BIT 8

// #define DEBUG

const float MAX_TRACE_DISTANCE = 15.0;
const float MIN_STEP_FACTOR = 0.2;
const float MIN_VIRTUAL_STEP_FACTOR = 0.01;
const float MIN_SPECULAR_APERTURE = 0.05;

uniform sampler2D u_diffuseTexture;
uniform sampler2D u_normalMap;
uniform sampler2D u_specularMap;
uniform sampler2D u_emissionMap;
uniform sampler2D u_depthTexture;
uniform sampler2D u_virtualMap;

uniform sampler3D u_voxelRadiance;
uniform sampler3D u_voxelOpacity;
uniform sampler3D u_voxelReflectance;
uniform sampler3D u_voxelNormal;
uniform sampler3D u_virtualVoxelRadiance;
uniform sampler3D u_virtualVoxelOpacity;
uniform sampler3D u_virtualVoxelDiffuse;
uniform sampler3D u_virtualVoxelSpecularA;
uniform sampler3D u_virtualVoxelNormal;

uniform float u_virtualVoxelSizeL0;
uniform vec3 u_virtualVolumeCenterL0;
uniform uint u_virtualVolumeDimension;

uniform int u_BRDFMode;
uniform mat4 u_viewProjInv;
uniform uint u_volumeDimension;
uniform vec3 u_eyePos;
uniform float u_voxelSizeL0;
uniform vec3 u_volumeCenterL0;
uniform float u_stepFactor;
uniform float u_viewAperture;
uniform int u_lightingMask;
uniform float u_indirectDiffuseIntensity;
uniform float u_indirectSpecularIntensity;
uniform float u_occlusionDecay = 1.0;
uniform float u_ambientOcclusionFactor;
uniform float u_traceStartOffset;
uniform int u_visualizeMinLevelSelection = 0;
uniform bool u_useNormalMapping;

// Debug
uniform float u_virtualStepFactor;
uniform int u_counterBreak;
uniform int u_virtualSelfOcclusion;
uniform float u_indirectSpecularShadow;
uniform float u_indirectDiffuseShadow;
uniform uint u_secondBounce;
uniform float u_secondIndirectDiffuse;
uniform uint u_realReflectance;

layout(location = 0) out vec4 out_color;

bool isVirtualFrag = uint(texture(u_virtualMap, In.texCoords).r) == 1;
// Variables dependent on isVirutal
vec3 frag_volumeCenterL0 = (isVirtualFrag)? u_virtualVolumeCenterL0 : u_volumeCenterL0;
float frag_voxelSizeL0 = (isVirtualFrag)? u_virtualVoxelSizeL0 : u_voxelSizeL0;
float frag_volumeDimension = (isVirtualFrag)? u_virtualVolumeDimension : u_volumeDimension;
uint frag_clipLevelCount = (isVirtualFrag)? VIRTUAL_CLIP_LEVEL_COUNT : CLIP_LEVEL_COUNT;
float frag_clipLevelCountInv = (isVirtualFrag)? VIRTUAL_CLIP_LEVEL_COUNT_INV : CLIP_LEVEL_COUNT_INV;

// #define USE_32_CONES

#ifdef USE_32_CONES
// 32 Cones for higher quality (16 on average per hemisphere)
const int DIFFUSE_CONE_COUNT = 32;
const float DIFFUSE_CONE_APERTURE = 0.628319; // 36 degree

const vec3 DIFFUSE_CONE_DIRECTIONS[32] = {
    vec3(0.898904, 0.435512, 0.0479745),
    vec3(0.898904, -0.435512, -0.0479745),
    vec3(0.898904, 0.0479745, -0.435512),
    vec3(0.898904, -0.0479745, 0.435512),
    vec3(-0.898904, 0.435512, -0.0479745),
    vec3(-0.898904, -0.435512, 0.0479745),
    vec3(-0.898904, 0.0479745, 0.435512),
    vec3(-0.898904, -0.0479745, -0.435512),
    vec3(0.0479745, 0.898904, 0.435512),
    vec3(-0.0479745, 0.898904, -0.435512),
    vec3(-0.435512, 0.898904, 0.0479745),
    vec3(0.435512, 0.898904, -0.0479745),
    vec3(-0.0479745, -0.898904, 0.435512),
    vec3(0.0479745, -0.898904, -0.435512),
    vec3(0.435512, -0.898904, 0.0479745),
    vec3(-0.435512, -0.898904, -0.0479745),
    vec3(0.435512, 0.0479745, 0.898904),
    vec3(-0.435512, -0.0479745, 0.898904),
    vec3(0.0479745, -0.435512, 0.898904),
    vec3(-0.0479745, 0.435512, 0.898904),
    vec3(0.435512, -0.0479745, -0.898904),
    vec3(-0.435512, 0.0479745, -0.898904),
    vec3(0.0479745, 0.435512, -0.898904),
    vec3(-0.0479745, -0.435512, -0.898904),
    vec3(0.57735, 0.57735, 0.57735),
    vec3(0.57735, 0.57735, -0.57735),
    vec3(0.57735, -0.57735, 0.57735),
    vec3(0.57735, -0.57735, -0.57735),
    vec3(-0.57735, 0.57735, 0.57735),
    vec3(-0.57735, 0.57735, -0.57735),
    vec3(-0.57735, -0.57735, 0.57735),
    vec3(-0.57735, -0.57735, -0.57735)
};
#else // 16 cones for lower quality (8 on average per hemisphere)
const int DIFFUSE_CONE_COUNT = 16;
const float DIFFUSE_CONE_APERTURE = 0.872665; // 50 degree

const vec3 DIFFUSE_CONE_DIRECTIONS[16] = {
    vec3(0.57735, 0.57735, 0.57735),
    vec3(0.57735, -0.57735, -0.57735),
    vec3(-0.57735, 0.57735, -0.57735),
    vec3(-0.57735, -0.57735, 0.57735),
    vec3(-0.903007, -0.182696, -0.388844),
    vec3(-0.903007, 0.182696, 0.388844),
    vec3(0.903007, -0.182696, 0.388844),
    vec3(0.903007, 0.182696, -0.388844),
    vec3(-0.388844, -0.903007, -0.182696),
    vec3(0.388844, -0.903007, 0.182696),
    vec3(0.388844, 0.903007, -0.182696),
    vec3(-0.388844, 0.903007, 0.182696),
    vec3(-0.182696, -0.388844, -0.903007),
    vec3(0.182696, 0.388844, -0.903007),
    vec3(-0.182696, 0.388844, 0.903007),
    vec3(0.182696, -0.388844, 0.903007)
};
#endif

vec3 worldPosFromDepth(float depth)
{
    vec4 p = vec4(In.texCoords, depth, 1.0);
    p.xyz = p.xyz * 2.0 - 1.0;
    p = u_viewProjInv * p;
    return p.xyz / p.w;
}

float getMinLevel(vec3 posW, vec3 volumeCenterL0, float voxelSizeL0, float volumeDimension)
{
    float distanceToCenter = length(volumeCenterL0 - posW);
    float minRadius = voxelSizeL0 * volumeDimension * 0.5;
    float minLevel = log2(distanceToCenter / minRadius);  
    minLevel = max(0.0, minLevel);

    float radius = minRadius * exp2(ceil(minLevel));
    float f = distanceToCenter / radius;
    
    // Smoothly transition from current level to the next level
    float transitionStart = 0.5;
    float c = 1.0 / (1.0 - transitionStart);
    
    return f > transitionStart ? ceil(minLevel) + (f - transitionStart) * c : ceil(minLevel);
}

bool inBox(vec3 p, AABBox3D box) {
    return all(greaterThan(p, box.minPos)) && all(lessThan(p, box.maxPos));
}

bool inSphere(vec3 p, vec3 center, float radius) {
    return length(p-center) < radius;
}

vec4 sampleVirtualClipmapTexture(sampler3D virtualClipmapTexture, vec3 posW, int clipmapLevel, vec3 faceOffsets, vec3 weight)
{
    float virtualExtent = u_virtualVoxelSizeL0 * exp2(clipmapLevel) * u_virtualVolumeDimension;

#ifdef VOXEL_TEXTURE_WITH_BORDER
	vec3 samplePos = (fract(posW / virtualExtent) * u_virtualVolumeDimension + vec3(BORDER_WIDTH)) / (float(u_virtualVolumeDimension) + 2.0 * BORDER_WIDTH);
#else
    vec3 samplePos = fract(posW / virtualExtent);
#endif

    samplePos.y += clipmapLevel;
    samplePos.y *= VIRTUAL_CLIP_LEVEL_COUNT_INV;
    samplePos.x *= FACE_COUNT_INV;
    return clamp(texture(virtualClipmapTexture, samplePos + vec3(faceOffsets.x, 0.0, 0.0)) * weight.x +
                 texture(virtualClipmapTexture, samplePos + vec3(faceOffsets.y, 0.0, 0.0)) * weight.y +
                 texture(virtualClipmapTexture, samplePos + vec3(faceOffsets.z, 0.0, 0.0)) * weight.z, 0.0, 1.0);
}

vec4 sampleVirtualClipmapLinearly(sampler3D clipmapTexture, vec3 posW, float curLevel, vec3 faceOffsets, vec3 weight)
{    
    int lowerLevel = int(floor(curLevel));
    int upperLevel = int(ceil(curLevel));
    
    vec4 lowSample = sampleVirtualClipmapTexture(clipmapTexture, posW, lowerLevel, faceOffsets, weight);
    
	if (lowerLevel == upperLevel)
        return lowSample;
	
    vec4 highSample = sampleVirtualClipmapTexture(clipmapTexture, posW, upperLevel, faceOffsets, weight);
	
    return mix(lowSample, highSample, fract(curLevel));
}

vec4 castConeVirtual(vec3 startPos, vec3 direction, float aperture, float maxDistance, float startLevel, 
                    out Intersection isect)
{
    isect.hit = false;
    float minRadius = u_virtualVoxelSizeL0 * u_virtualVolumeDimension * 0.5;

    float curLevel = startLevel;
    float voxelSize = u_virtualVoxelSizeL0 * exp2(curLevel);

    // 2 is multiplied to avoid missing voxel region for higher clipmap level
    float maxVolumeExtent = u_virtualVoxelSizeL0 * u_virtualVolumeDimension * 2;

    // Find the intersection point with the AABB of virtual objectfor at the cone direction
    // and use the point as a start startLevel for cone casting
    float s = 0.0;
    Ray ray = Ray(startPos, direction);
    AABBox3D aabb = AABBox3D(u_virtualVolumeCenterL0 - maxVolumeExtent * 0.5f, u_virtualVolumeCenterL0 + maxVolumeExtent * 0.5f);
    vec2 tMinMax = rayIntersectsAABB(ray, aabb);
    bool intersected = (tMinMax[0] < tMinMax[1] && tMinMax[0] > 0.0) || inBox(startPos, aabb);
    if (!intersected) return vec4(vec3(0), 1); // It means no occlusion

    // Start from the object's bounding box if the startPos is not in the box
    // and end at the object's bounding box
    if (tMinMax[0] > 0.f) {
        s = tMinMax[0];
    }
    maxDistance = tMinMax[1];

    // Initialize accumulated color and opacity
    vec4 dst = vec4(0.0);
    // Coefficient used in the computation of the diameter of a cone
	float coneCoefficient = 2.0 * tan(aperture * 0.5);

    // Diameter of cone at start position is the l0 voxel size
    float diameter = max(s * coneCoefficient, u_virtualVoxelSizeL0);

    float stepFactor = max(MIN_VIRTUAL_STEP_FACTOR, u_virtualStepFactor);
	float occlusion = 0.0;
    
    ivec3 faceIndices = computeVoxelFaceIndices(direction); // Implementation in voxelConeTracing/common.glsl
    vec3 faceOffsets = vec3(faceIndices) * FACE_COUNT_INV;
    vec3 weight = direction * direction;
    
    float curSegmentLength = voxelSize;

    int counter = 0;
    // Ray marching - compute occlusion and radiance in one go
    while (s < maxDistance && occlusion < 1.0)
    {
        vec3 position = startPos + direction * s;

        float distanceToCenter = length(u_virtualVolumeCenterL0 - position);
        float minLevel = ceil(log2(distanceToCenter / minRadius));
        
        curLevel = log2(diameter / u_virtualVoxelSizeL0);
        curLevel = min(max(max(startLevel, curLevel), minLevel), VIRTUAL_CLIP_LEVEL_COUNT - 1);
        
        // Retrieve radiance by accessing the 3D clipmap (voxel radiance and opacity)
        vec4 radiance;
        
        // Fragment for real object
        //  - Use virtual map for shadow (differential)
        vec4 virtualOpacity = sampleVirtualClipmapLinearly(u_virtualVoxelOpacity, position, curLevel, faceOffsets, weight);
        radiance = vec4(vec3(0), virtualOpacity.a); // virtualOpacity.rgb is for second bounce (material, normal, etc.)
		float opacity = radiance.a;
        if (!isect.hit && opacity > EPSILON) {
            isect.position = position;
            isect.normal = unpackNormal(sampleVirtualClipmapLinearly(u_virtualVoxelNormal, position, curLevel, faceIndices, weight).rgb);
            isect.level = curLevel;
            isect.hit = true;
            if(isect.hitOnly) break;
        }

        voxelSize = u_virtualVoxelSizeL0 * exp2(curLevel);
        
        // Radiance correction
        float correctionQuotient = curSegmentLength / voxelSize;
        radiance.rgb = radiance.rgb * correctionQuotient;
		
        // Opacity correction
        opacity = clamp(1.0 - pow(1.0 - opacity, correctionQuotient), 0.0, 1.0);

        vec4 src = vec4(radiance.rgb, opacity);
		
        // Front-to-back compositing
        dst += clamp(1.0 - dst.a, 0.0, 1.0) * src;
		occlusion += (1.0 - occlusion) * opacity / (1.0 + (s + voxelSize) * u_occlusionDecay);

		float sLast = s;
        s += max(diameter, u_virtualVoxelSizeL0) * stepFactor;
        curSegmentLength = (s - sLast);
        diameter = s * coneCoefficient;

        if (counter++ > u_counterBreak) break;
    }

    return clamp(vec4(dst.rgb, 1.0 - occlusion), 0.0, 1.0);
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
    
    vec4 lowSample = sampleClipmapTexture(clipmapTexture, posW, lowerLevel, faceOffsets, weight);
    
	if (lowerLevel == upperLevel)
        return lowSample;
	
    vec4 highSample = sampleClipmapTexture(clipmapTexture, posW, upperLevel, faceOffsets, weight);
	
    return mix(lowSample, highSample, fract(curLevel));
}

vec4 castCone(vec3 startPos, vec3 direction, float aperture, float maxDistance, float startLevel, 
            out Intersection isect)
{
    isect.hit = false;
    // Initialize accumulated color and opacity
    vec4 dst = vec4(0.0);
    // Coefficient used in the computation of the diameter of a cone
	float coneCoefficient = 2.0 * tan(aperture * 0.5);
    
    float curLevel = startLevel;
    float voxelSize = u_voxelSizeL0 * exp2(curLevel);
    
    // Offset startPos in the direction to avoid self occlusion and reduce voxel aliasing
    startPos += direction * voxelSize * u_traceStartOffset * 0.5;

    // Distance from ray origin to current step
    float s = 0.0;
    
    // Diameter of cone at s
    float diameter = max(s * coneCoefficient, u_voxelSizeL0);

    float stepFactor = max(MIN_STEP_FACTOR, u_stepFactor);
	float occlusion = 0.0;
    
    ivec3 faceIndices = computeVoxelFaceIndices(direction); // Implementation in voxelConeTracing/common.glsl
    vec3 weight = direction * direction;
    
    float curSegmentLength = voxelSize;
    
    // minimum radius from clip region center to position
    float minRadius = u_voxelSizeL0 * u_volumeDimension * 0.5;
    
    // Ray marching - compute occlusion and radiance in one go
    while (s < maxDistance && occlusion < 1.0)
    {
        vec3 position = startPos + direction * s;
        
        float distanceToCenter = length(u_volumeCenterL0 - position);
        float minLevel = ceil(log2(distanceToCenter / minRadius));
        
        curLevel = log2(diameter / u_voxelSizeL0);
        // The startLevel is the minimum level we start off with, minLevel is the current minLevel
        // It's important to use the max of both (and curLevel of course) because we don't want to suddenly
        // sample at a lower level than we started off with and ensure that we don't sample in a level that is too low.
        curLevel = min(max(max(startLevel, curLevel), minLevel), CLIP_LEVEL_COUNT - 1);
        
        // Retrieve radiance by accessing the 3D clipmap (voxel radiance and opacity)
        vec4 radiance = sampleClipmapLinearly(u_voxelRadiance, position, curLevel, faceIndices, weight);
		float opacity = radiance.a;
        if (!isect.hit && opacity > EPSILON) {
            isect.position = position;
            isect.normal = unpackNormal(sampleClipmapLinearly(u_voxelNormal, position, curLevel, faceIndices, weight).rgb);
            isect.level = curLevel;
            isect.hit = true;
            if(isect.hitOnly) break;
        }

        voxelSize = u_voxelSizeL0 * exp2(curLevel);
        
        // Radiance correction
        float correctionQuotient = curSegmentLength / voxelSize;
        radiance.rgb = radiance.rgb * correctionQuotient;
		
        // Opacity correction
        opacity = clamp(1.0 - pow(1.0 - opacity, correctionQuotient), 0.0, 1.0);

        vec4 src = vec4(radiance.rgb, opacity);
		
        // Front-to-back compositing
        dst += clamp(1.0 - dst.a, 0.0, 1.0) * src;
		occlusion += (1.0 - occlusion) * opacity / (1.0 + (s + voxelSize) * u_occlusionDecay);

		float sLast = s;
        s += max(diameter, u_voxelSizeL0) * stepFactor;
        curSegmentLength = (s - sLast);
        diameter = s * coneCoefficient;
    }
    
    return clamp(vec4(dst.rgb, 1.0 - occlusion), 0.0, 1.0);
}

// Used to visualize the transition of minLevel selection
vec4 minLevelToColor(float minLevel)
{
   vec4 colors[] = {
        vec4(1.0, 0.0, 0.0, 1.0),
        vec4(0.0, 1.0, 0.0, 1.0),
        vec4(0.0, 0.0, 1.0, 1.0),
        vec4(1.0, 1.0, 0.0, 1.0),
        vec4(0.0, 1.0, 1.0, 1.0),
        vec4(1.0, 0.0, 1.0, 1.0),
        vec4(1.0, 1.0, 1.0, 1.0)
    };
	
	vec4 minLevelColor = vec4(0.0);
    
    if (minLevel < 1)
        minLevelColor = mix(colors[0], colors[1], fract(minLevel));
    else if (minLevel < 2)
        minLevelColor = mix(colors[1], colors[2], fract(minLevel));
    else if (minLevel < 3)
        minLevelColor = mix(colors[2], colors[3], fract(minLevel));
    else if (minLevel < 4)
        minLevelColor = mix(colors[3], colors[4], fract(minLevel));
    else if (minLevel < 5)
        minLevelColor = mix(colors[4], colors[5], fract(minLevel));
    else if (minLevel < 6)
        minLevelColor = mix(colors[5], colors[6], fract(minLevel));
        
	return minLevelColor * 0.5;
}

vec4 castDiffuseCones(vec3 startPos, vec3 normal, float minLevel, bool traceVirtual, bool inverse = false)
{
    Intersection tmpIsect;
    vec4 indirectContribution = vec4(0.0);
    for (int i = 0; i < DIFFUSE_CONE_COUNT; ++i)
    {
		float cosTheta = dot(normal, DIFFUSE_CONE_DIRECTIONS[i]);
        
        if (cosTheta < 0.0)
            continue;

        if (traceVirtual && inverse) {
            indirectContribution += (1 - castConeVirtual(startPos, DIFFUSE_CONE_DIRECTIONS[i], DIFFUSE_CONE_APERTURE, MAX_TRACE_DISTANCE, minLevel, tmpIsect)) * cosTheta;
        // 1nd: from real to virtual
        } else if (traceVirtual && !inverse) {

            // Find hitpoint with virtual object
            Intersection isect;
            indirectContribution += castConeVirtual(startPos, DIFFUSE_CONE_DIRECTIONS[i], DIFFUSE_CONE_APERTURE, MAX_TRACE_DISTANCE, minLevel, isect) * cosTheta;

            // vec3 weight = DIFFUSE_CONE_DIRECTIONS[i] * DIFFUSE_CONE_DIRECTIONS[i];
            // ivec3 outFaceIndices = computeVoxelFaceIndices(-DIFFUSE_CONE_DIRECTIONS[i]);
            // vec3 reflectance = sampleClipmapLinearly(u_voxelReflectance, startPos, minLevel, outFaceIndices, weight).rgb;
            // indirectContribution.rgb += reflectance;
            // continue;

            // if (isect.hit) {
            //     indirectContribution.rgb = isect.position;
            //     return indirectContribution;
            // }

            // 2nd: from virtual to real
            if (isect.hit && u_secondBounce == 1) {
                vec3 weight = DIFFUSE_CONE_DIRECTIONS[i] * DIFFUSE_CONE_DIRECTIONS[i];
                ivec3 inFaceIndices = computeVoxelFaceIndices(DIFFUSE_CONE_DIRECTIONS[i]);
                // minLevel: virtualMinLevel
                vec3 diffuse = sampleVirtualClipmapLinearly(u_virtualVoxelDiffuse, isect.position, minLevel, inFaceIndices, weight).rgb;
                vec4 specularA = sampleVirtualClipmapLinearly(u_virtualVoxelSpecularA, isect.position, minLevel, inFaceIndices, weight);
                vec3 specular = specularA.rgb;
                float shininess = specularA.a;
                float roughness = shininessToRoughness(shininess);

                // minLevel2: real min level at hit position
                float minLevel2 = getMinLevel(isect.position, u_volumeCenterL0, u_voxelSizeL0, u_volumeDimension);   
                vec4 diffuseContribution = vec4(0.0);
                for (int j = 0; j < DIFFUSE_CONE_COUNT; ++j) {
                    float cosTheta2 = dot(isect.normal, DIFFUSE_CONE_DIRECTIONS[j]);
                    
                    if (cosTheta2 < 0.0)
                        continue;

                    diffuseContribution += castCone(isect.position, DIFFUSE_CONE_DIRECTIONS[j], DIFFUSE_CONE_APERTURE, MAX_TRACE_DISTANCE, minLevel2, tmpIsect) * cosTheta2;
                }
                diffuseContribution /= DIFFUSE_CONE_COUNT * 0.5;
                indirectContribution.rgb += diffuse * diffuseContribution.rgb * u_secondIndirectDiffuse;
            }
        // : from virtual to real
        } else { 
            indirectContribution += castCone(startPos, DIFFUSE_CONE_DIRECTIONS[i], DIFFUSE_CONE_APERTURE, MAX_TRACE_DISTANCE, minLevel, tmpIsect) * cosTheta;
        }
    }
    return indirectContribution;
}

void main() 
{
    float depth = texture2D(u_depthTexture, In.texCoords).r;
    if (depth == 1.0)
        discard;
    // if (isVirtualFrag) discard;
    vec3 diffuse = texture(u_diffuseTexture, In.texCoords).rgb;
    vec3 normal = unpackNormal(texture(u_normalMap, In.texCoords).rgb);
    vec3 emission = texture(u_emissionMap, In.texCoords).rgb;
    bool hasEmission = any(greaterThan(emission, vec3(0.0)));  
    vec3 posW = worldPosFromDepth(depth);
    vec3 view = normalize(u_eyePos - posW);
    vec4 specColor = texture(u_specularMap, In.texCoords);
    
    float minLevel = getMinLevel(posW, u_volumeCenterL0, u_voxelSizeL0, u_volumeDimension);
    float virtualMinLevel = getMinLevel(posW, u_virtualVolumeCenterL0, u_virtualVoxelSizeL0, u_virtualVolumeDimension);
    
    // Compute indirect contribution
    vec4 indirectContribution = vec4(0.0);
    
    // Offset the starting pos to avoid the self-occlusion
    float voxelSize = frag_voxelSizeL0 * exp2(isVirtualFrag? virtualMinLevel : minLevel);
    vec3 startPosOffset = posW + normal * voxelSize * u_traceStartOffset;
    vec3 startPos = posW;

	if (isVirtualFrag){
        // Radiance
        indirectContribution = castDiffuseCones(startPos, normal, minLevel, false);
        // Delta (shadow)
        if (u_virtualSelfOcclusion == 1) {
            float a = castDiffuseCones(startPosOffset, normal, virtualMinLevel, true, true).a;
            indirectContribution.rgb -= vec3(a) * u_indirectDiffuseShadow;
        }
    }
    else {
        ivec3 outFaceIndices = computeVoxelFaceIndices(-normal);
        vec3 weight = normal * normal;
        vec3 reflectance = sampleClipmapLinearly(u_voxelReflectance, startPos, minLevel, outFaceIndices, weight).rgb;
        // out_color = vec4(reflectance, 1);
        // return;

        // Occlusion-based shadow + color bleeding of virtual object
        indirectContribution = castDiffuseCones(startPos, normal, virtualMinLevel, true);
        indirectContribution.rgb *= reflectance;
    }

    // DIFFUSE_CONE_COUNT includes cones to integrate over a sphere - on the hemisphere there are on average ~half of these cones
	indirectContribution /= DIFFUSE_CONE_COUNT * 0.5;
    indirectContribution.a *= u_ambientOcclusionFactor;
#ifdef DEBUG
    out_color = vec4(vec3(indirectContribution.rgb), 1) * u_indirectDiffuseIntensity;
    return;
#endif
    
	indirectContribution.rgb *= diffuse * u_indirectDiffuseIntensity;
    indirectContribution = clamp(indirectContribution, 0.0, 1.0);
    
	// Specular cone
    vec3 specularConeDirection = reflect(-view, normal);
    vec3 specularContribution = vec3(0.0);
    
    float shininess = unpackShininess(specColor.a);
    float roughness = shininessToRoughness(shininess);
    
    // Currentlly, only for the (isVirtualFrag == true)
    if (any(greaterThan(specColor.rgb, vec3(EPSILON))) && specColor.a > EPSILON) {
        // Radiance
        Intersection isect;
        specularContribution = castCone(startPos, specularConeDirection, max(roughness, MIN_SPECULAR_APERTURE), MAX_TRACE_DISTANCE, minLevel, isect).rgb * specColor.rgb * u_indirectSpecularIntensity;

        // Delta
        if (u_virtualSelfOcclusion == 1) {
            float a = castConeVirtual(startPosOffset, specularConeDirection, max(roughness, MIN_SPECULAR_APERTURE), MAX_TRACE_DISTANCE, virtualMinLevel, isect).a;
            // anti radiance
            specularContribution -= vec3(1-a) * specColor.rgb * u_indirectSpecularShadow;
        }
        
        // Find hit point for virtual object and cast difffuse cones again
        
    }
    
    vec3 directContribution = vec3(0.0);
    
    if (hasEmission)
    {
        directContribution += emission;
    }
    
    directContribution = clamp(directContribution, 0.0, 1.0);
    
	out_color = vec4(0.0, 0.0, 0.0, 1.0);
    
    if ((u_lightingMask & AMBIENT_OCCLUSION_BIT) != 0)
        directContribution *= indirectContribution.a;

    if ((u_lightingMask & DIRECT_LIGHTING_BIT) != 0)
        out_color.rgb += directContribution;
        
    if ((u_lightingMask & INDIRECT_DIFFUSE_LIGHTING_BIT) != 0)
        out_color.rgb += indirectContribution.rgb;
        
    if ((u_lightingMask & INDIRECT_SPECULAR_LIGHTING_BIT) != 0)
        out_color.rgb += specularContribution;
        
    // If only ambient occlusion is selected show ambient occlusion
    if (u_lightingMask == AMBIENT_OCCLUSION_BIT)
        out_color.rgb = vec3(indirectContribution.a);
        
    if (u_visualizeMinLevelSelection > 0)
        out_color *= minLevelToColor(minLevel);

    if (isVirtualFrag) {
        out_color = clamp(out_color, 0.0, 1.0);
    }
    else {
        // original shading (VCT)
        out_color = clamp(out_color, 0.0, 1.0);
    }
}
