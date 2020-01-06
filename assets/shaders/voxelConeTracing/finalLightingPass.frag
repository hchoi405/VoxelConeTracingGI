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

#define DIRECT_LIGHTING_BIT 1
#define INDIRECT_DIFFUSE_LIGHTING_BIT 2
#define INDIRECT_SPECULAR_LIGHTING_BIT 4
#define AMBIENT_OCCLUSION_BIT 8

// #define DEBUG

const float MAX_TRACE_DISTANCE = 30.0;
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
uniform sampler3D u_virtualVoxelRadiance;
uniform sampler3D u_virtualVoxelOpacity;

uniform float u_virtualVoxelSizeL0;
uniform vec3 u_virtualVolumeCenterL0;
uniform uint u_virtualVolumeDimension;

uniform DirectionalLight u_directionalLights[MAX_DIR_LIGHT_COUNT];
uniform DirectionalLightShadowDesc u_directionalLightShadowDescs[MAX_DIR_LIGHT_COUNT];
uniform sampler2D u_shadowMaps[MAX_DIR_LIGHT_COUNT];
uniform int u_numActiveDirLights;
uniform float u_depthBias;
uniform float u_usePoissonFilter;

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

bool inBox(vec3 p, vec3 low, vec3 high) {
    return all(greaterThan(p, low)) && all(lessThan(p, high));
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

vec4 castConeVirtual(vec3 startPos, vec3 direction, float aperture, float maxDistance, float startLevel)
{
    // float distanceToCenter = length(u_virtualVolumeCenterL0 - startPos);
    float minRadius = u_virtualVoxelSizeL0 * u_virtualVolumeDimension * 0.5;
    // float minLevel = log2(distanceToCenter / minRadius);  
    // startLevel = max(0.0, minLevel);

    // Find the intersection point with the AABB of virtual objectfor at the cone direction
    // and use the point as a start startLevel for cone casting
    // Ray ray = Ray(startPos, direction);
    // AABBox3D aabb = AABBox3D(u_virtualMin, u_virtualMax);
    // vec2 tMinMax = rayIntersectsAABB(ray, aabb);
    // bool intersected = (tMinMax[0] < tMinMax[1] && tMinMax[0] > 0.0);
    bool intersected = false;

    // bool intersected = (tMinMax[0] < tMinMax[1] && tMinMax[0] > 0.0);
    // if (tMinMax[0] > tMinMax[1]) return vec4(0, 0, 1, 1);
    // if (!intersected) return vec4(0, 0, 0, 0);
    // if (tMinMax[0] < 0.0) return vec4(1, 0, 0, 1);
    // else return vec4(vec3(1),1);
    // else return vec4(1-vec3(tMinMax[0]) / u_indirectDiffuseIntensity, 1);

    // Initialize accumulated color and opacity
    vec4 dst = vec4(0.0);
    // Coefficient used in the computation of the diameter of a cone
	float coneCoefficient = 2.0 * tan(aperture * 0.5);
    
    // Start from the object's bounding box and end at the object's bounding box
    float s = 0.0;
    // float s = (intersected)? tMinMax[0] : 0.0;
    // maxDistance = (intersected)? tMinMax[1] : maxDistance;
    // return vec4(vec3(s), 1.0);

    // Diameter of cone at start position is the l0 voxel size
    float diameter = max(s * coneCoefficient, u_virtualVoxelSizeL0);

    // Skip too far regions
    // if (startLevel > VIRTUAL_CLIP_LEVEL_COUNT) return vec4(1,0,0,1);
    // else if (startLevel == VIRTUAL_CLIP_LEVEL_COUNT) return vec4(0,1,0,1);
    // else return vec4(0,0,1,1);
    // return vec4(vec3((VIRTUAL_CLIP_LEVEL_COUNT-startLevel) * VIRTUAL_CLIP_LEVEL_COUNT_INV),1);

    float curLevel = startLevel;
    float voxelSize = u_virtualVoxelSizeL0 * exp2(curLevel);

    float stepFactor = max(MIN_VIRTUAL_STEP_FACTOR, u_virtualStepFactor);
	float occlusion = 0.0;
    
    ivec3 faceIndices = computeVoxelFaceIndices(direction); // Implementation in voxelConeTracing/common.glsl
    vec3 faceOffsets = vec3(faceIndices) * FACE_COUNT_INV;
    vec3 weight = direction * direction;
    
    float curSegmentLength = voxelSize;

    // float minRadius = u_virtualVoxelSizeL0 * u_virtualVolumeDimension * 0.5;

    // State variables to check wether hit virtual or real object
    bool hitVirtual = false;
    bool hitReal = false;

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

vec4 castCone(vec3 startPos, vec3 direction, float aperture, float maxDistance, float startLevel)
{
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

void main() 
{
    float depth = texture2D(u_depthTexture, In.texCoords).r;
    if (depth == 1.0)
        discard;
        
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
    
    // TODO: Check if offset is needed
    float voxelSize = frag_voxelSizeL0 * exp2(isVirtualFrag? virtualMinLevel : minLevel);
    vec3 startPosOffset = posW + normal * voxelSize * u_traceStartOffset;
    vec3 startPos = posW;
    
    float cosSum = 0.0;
	for (int i = 0; i < DIFFUSE_CONE_COUNT; ++i)
    {
		float cosTheta = dot(normal, DIFFUSE_CONE_DIRECTIONS[i]);
        
        if (cosTheta < 0.0)
            continue;
        
        if (isVirtualFrag){
            float a = castConeVirtual(startPosOffset, DIFFUSE_CONE_DIRECTIONS[i], DIFFUSE_CONE_APERTURE, MAX_TRACE_DISTANCE, virtualMinLevel).a;
            indirectContribution += castCone(startPos, DIFFUSE_CONE_DIRECTIONS[i], DIFFUSE_CONE_APERTURE ,MAX_TRACE_DISTANCE, minLevel) * cosTheta;
            indirectContribution.rgb -= vec3(1-a) * u_indirectDiffuseShadow * cosTheta;
        }
        else {
            indirectContribution += castConeVirtual(startPos, DIFFUSE_CONE_DIRECTIONS[i], DIFFUSE_CONE_APERTURE, MAX_TRACE_DISTANCE, virtualMinLevel) * cosTheta;
        }
    }

    // DIFFUSE_CONE_COUNT includes cones to integrate over a sphere - on the hemisphere there are on average ~half of these cones
	indirectContribution /= DIFFUSE_CONE_COUNT * 0.5;
    indirectContribution.a *= u_ambientOcclusionFactor;
#ifdef DEBUG
    out_color = vec4(vec3(indirectContribution.a), 1);
    return;
#endif
    
	indirectContribution.rgb *= diffuse * u_indirectDiffuseIntensity;
    indirectContribution = clamp(indirectContribution, 0.0, 1.0);
    
	// Specular cone
    vec3 specularConeDirection = reflect(-view, normal);
    vec3 specularContribution = vec3(0.0);
    
    float shininess = unpackShininess(specColor.a);
    float roughness = shininessToRoughness(shininess);
    
    if (any(greaterThan(specColor.rgb, vec3(EPSILON))) && specColor.a > EPSILON) {
        specularContribution = castCone(startPos, specularConeDirection, max(roughness, MIN_SPECULAR_APERTURE), MAX_TRACE_DISTANCE, minLevel).rgb * specColor.rgb * u_indirectSpecularIntensity;
        float a = castConeVirtual(startPosOffset, specularConeDirection, max(roughness, MIN_SPECULAR_APERTURE), MAX_TRACE_DISTANCE, virtualMinLevel).a;
        if (u_virtualSelfOcclusion == 1)
            specularContribution -= vec3(1-a) * specColor.rgb * u_indirectSpecularShadow;
    }
    
    vec3 directContribution = vec3(0.0);
    
    if (hasEmission)
    {
        directContribution += emission;
    }
    else
    {
        for (int i = 0; i < u_numActiveDirLights; ++i)
        {
            vec3 lightDir = u_directionalLights[i].direction;
            float nDotL = max(0.0, dot(normal, -lightDir));
            vec3 halfway = normalize(view - lightDir);
            
            float visibility = 1.0;
            if (u_directionalLightShadowDescs[i].enabled != 0)
            {
                visibility = computeVisibility(posW, u_shadowMaps[i], u_directionalLightShadowDescs[i], u_usePoissonFilter, u_depthBias);
            }
            
            vec3 lightColor = u_directionalLights[i].color;
            
            if (u_BRDFMode == BLINN_PHONG_MODE_IDX)
            {
                vec3 blinnPhong = blinnPhongBRDF(lightColor, diffuse, lightColor, specColor.rgb, normal, -lightDir, halfway, shininess);
                directContribution += visibility * blinnPhong * u_directionalLights[i].intensity;
            } else if (u_BRDFMode == COOK_TORRANCE_MODE_IDX)
            {
                vec3 cook = cookTorranceBRDF(-lightDir, normal, view, halfway, roughness, specColor.rgb * 0.5);
                directContribution += visibility * (cook * lightColor * specColor.rgb + lightColor * diffuse * nDotL) * u_directionalLights[i].intensity;
            }
        }
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

        // max aperture
        // out_color = vec4(castCone(u_eyePos, -view, sqrt(2/3), MAX_TRACE_DISTANCE, getMinLevel(u_eyePos)).rgb, 1.0) * u_indirectSpecularIntensity;

        // control aperture
        // out_color = vec4(castCone(u_eyePos, -view, u_viewAperture, MAX_TRACE_DISTANCE, 0).rgb, 1.0) * u_indirectSpecularIntensity;

        // direct evaluation
        // ivec3 faceIndices = computeVoxelFaceIndices(-view);
        // vec3 w = view * view;
        // vec3 w = vec3(1.0);
        // out_color = sampleClipmapLinearly(u_voxelRadiance, posW, minLevel, faceIndices, w);
    }
}
