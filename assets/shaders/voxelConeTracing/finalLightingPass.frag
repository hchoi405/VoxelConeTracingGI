#version 430
#extension GL_ARB_shading_language_include : enable

#include "/voxelConeTracing/settings.glsl"
#include "/voxelConeTracing/common.glsl"
#include "/BRDF.glsl"
#include "/shadows/shadows.glsl"
#include "/intersection.glsl"
#include "/voxelConeTracing/conversion.glsl"
#include "/commonStruct.h"
#include "/util/noise.glsl"

in Vertex { vec2 texCoords; }
In;

#define DIRECT_LIGHTING_BIT 1
#define INDIRECT_DIFFUSE_LIGHTING_BIT 2
#define INDIRECT_SPECULAR_LIGHTING_BIT 4
#define AMBIENT_OCCLUSION_BIT 8

const float MAX_TRACE_DISTANCE = 15.0;
const float MIN_STEP_FACTOR = 0.2;
// const float MIN_SPECULAR_APERTURE = 0.1; // 5.73 degrees
const float MIN_SPECULAR_APERTURE = 0.05;  // 2.86 degrees

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
uniform float u_virtualStepFactor;
uniform float u_viewAperture;
uniform int u_lightingMask;
uniform float u_realIndirectDiffuseIntensity;
uniform float u_virtualIndirectDiffuseIntensity;
uniform float u_indirectSpecularIntensity;
uniform float u_occlusionDecay = 1.0;
uniform float u_ambientOcclusionFactor;
uniform float u_traceStartOffset;
uniform float u_traceDirectionOffset;
uniform int u_visualizeMinLevelSelection = 0;
uniform bool u_useNormalMapping;

sampler3D getOpacity(bool isVirtual) { return isVirtual ? u_virtualVoxelOpacity : u_voxelOpacity; }
sampler3D getRadiance(bool isVirtual) { return isVirtual ? u_virtualVoxelDiffuse : u_voxelRadiance; }
sampler3D getNormal(bool isVirtual) { return isVirtual ? u_virtualVoxelNormal : u_voxelNormal; }
sampler3D getDiffuse(bool isVirtual) { return isVirtual ? u_virtualVoxelDiffuse : u_voxelReflectance; }
sampler3D getSpecularA(bool isVirtual) { return isVirtual ? u_virtualVoxelSpecularA : u_voxelReflectance; }

VCTScene virtualScene, realScene;
VCTParams params;

// Debug
uniform uint u_debugFlag;
uniform uint u_renderReal;
uniform uint u_renderVirtual;
uniform int u_virtualSelfOcclusion;
uniform float u_indirectSpecularShadow;
uniform float u_indirectDiffuseShadow;
uniform uint u_secondBounce;
uniform float u_secondIndirectDiffuseIntensity;
uniform float u_secondIndirectSpecularIntensity;
uniform uint u_realReflectance;
uniform float u_ambientSecondIntensity;

layout(location = 0) out vec4 out_color;

bool isVirtualFrag = uint(texture(u_virtualMap, In.texCoords).r) == 1;
// Variables dependent on isVirutal
float frag_voxelSizeL0 = (isVirtualFrag) ? u_virtualVoxelSizeL0 : u_voxelSizeL0;

// #define USE_32_CONES

#ifdef USE_32_CONES
// 32 Cones for higher quality (16 on average per hemisphere)
const int DIFFUSE_CONE_COUNT = 32;
const float DIFFUSE_CONE_APERTURE = 0.628319;  // 36 degree

const vec3 DIFFUSE_CONE_DIRECTIONS[32] = {vec3(0.898904, 0.435512, 0.0479745),   vec3(0.898904, -0.435512, -0.0479745),
                                          vec3(0.898904, 0.0479745, -0.435512),  vec3(0.898904, -0.0479745, 0.435512),
                                          vec3(-0.898904, 0.435512, -0.0479745), vec3(-0.898904, -0.435512, 0.0479745),
                                          vec3(-0.898904, 0.0479745, 0.435512),  vec3(-0.898904, -0.0479745, -0.435512),
                                          vec3(0.0479745, 0.898904, 0.435512),   vec3(-0.0479745, 0.898904, -0.435512),
                                          vec3(-0.435512, 0.898904, 0.0479745),  vec3(0.435512, 0.898904, -0.0479745),
                                          vec3(-0.0479745, -0.898904, 0.435512), vec3(0.0479745, -0.898904, -0.435512),
                                          vec3(0.435512, -0.898904, 0.0479745),  vec3(-0.435512, -0.898904, -0.0479745),
                                          vec3(0.435512, 0.0479745, 0.898904),   vec3(-0.435512, -0.0479745, 0.898904),
                                          vec3(0.0479745, -0.435512, 0.898904),  vec3(-0.0479745, 0.435512, 0.898904),
                                          vec3(0.435512, -0.0479745, -0.898904), vec3(-0.435512, 0.0479745, -0.898904),
                                          vec3(0.0479745, 0.435512, -0.898904),  vec3(-0.0479745, -0.435512, -0.898904),
                                          vec3(0.57735, 0.57735, 0.57735),       vec3(0.57735, 0.57735, -0.57735),
                                          vec3(0.57735, -0.57735, 0.57735),      vec3(0.57735, -0.57735, -0.57735),
                                          vec3(-0.57735, 0.57735, 0.57735),      vec3(-0.57735, 0.57735, -0.57735),
                                          vec3(-0.57735, -0.57735, 0.57735),     vec3(-0.57735, -0.57735, -0.57735)};
#else  // 16 cones for lower quality (8 on average per hemisphere)
const int DIFFUSE_CONE_COUNT = 16;
const float DIFFUSE_CONE_APERTURE = 0.872665;  // 50 degree

const vec3 DIFFUSE_CONE_DIRECTIONS[16] = {
    vec3(0.57735, 0.57735, 0.57735),       vec3(0.57735, -0.57735, -0.57735),     vec3(-0.57735, 0.57735, -0.57735),
    vec3(-0.57735, -0.57735, 0.57735),     vec3(-0.903007, -0.182696, -0.388844), vec3(-0.903007, 0.182696, 0.388844),
    vec3(0.903007, -0.182696, 0.388844),   vec3(0.903007, 0.182696, -0.388844),   vec3(-0.388844, -0.903007, -0.182696),
    vec3(0.388844, -0.903007, 0.182696),   vec3(0.388844, 0.903007, -0.182696),   vec3(-0.388844, 0.903007, 0.182696),
    vec3(-0.182696, -0.388844, -0.903007), vec3(0.182696, 0.388844, -0.903007),   vec3(-0.182696, 0.388844, 0.903007),
    vec3(0.182696, -0.388844, 0.903007)};
#endif

vec3 worldPosFromDepth(float depth) {
    vec4 p = vec4(In.texCoords, depth, 1.0);
    p.xyz = p.xyz * 2.0 - 1.0;
    p = u_viewProjInv * p;
    return p.xyz / p.w;
}

float getMinLevel(vec3 posW, vec3 volumeCenterL0, float voxelSizeL0, float volumeDimension) {
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

float getRealMinLevel(vec3 posW) { return getMinLevel(posW, u_volumeCenterL0, u_voxelSizeL0, u_volumeDimension); }

float getVirtualMinLevel(vec3 posW) {
    return getMinLevel(posW, u_virtualVolumeCenterL0, u_virtualVoxelSizeL0, u_virtualVolumeDimension);
}

bool inBox(vec3 p, AABBox3D box) { return all(greaterThan(p, box.minPos)) && all(lessThan(p, box.maxPos)); }

bool inSphere(vec3 p, vec3 center, float radius) { return length(p - center) < radius; }

vec4 sampleUnifiedClipmapTexture(sampler3D virtualClipmapTexture, vec3 posW, int clipmapLevel, float voxelSizeL0,
                                 uint volumeDimension, float maxCliplevelInv, vec3 faceOffsets, vec3 weight) {
    float voxelSize = voxelSizeL0 * exp2(clipmapLevel);
    float extent = voxelSize * volumeDimension;

#ifdef VOXEL_TEXTURE_WITH_BORDER
    vec3 samplePos =
        (fract(posW / extent) * volumeDimension + vec3(BORDER_WIDTH)) / (float(volumeDimension) + 2.0 * BORDER_WIDTH);
#else
    vec3 samplePos = fract(posW / extent);
#endif

    samplePos.y += clipmapLevel;
    samplePos.y *= maxCliplevelInv;
    samplePos.x *= FACE_COUNT_INV;
    return clamp(texture(virtualClipmapTexture, samplePos + vec3(faceOffsets.x, 0.0, 0.0)) * weight.x +
                     texture(virtualClipmapTexture, samplePos + vec3(faceOffsets.y, 0.0, 0.0)) * weight.y +
                     texture(virtualClipmapTexture, samplePos + vec3(faceOffsets.z, 0.0, 0.0)) * weight.z,
                 0.0, 1.0);
}

vec4 sampleUnifiedClipmapLinearly(sampler3D clipmapTexture, vec3 position, float curLevel, float voxelSizeL0,
                                  uint volumeDimension, float maxCliplevelInv, vec3 faceOffsets, vec3 weight) {
    int lowerLevel = int(floor(curLevel));
    int upperLevel = int(ceil(curLevel));

    vec4 lowSample = sampleUnifiedClipmapTexture(clipmapTexture, position, lowerLevel, voxelSizeL0, volumeDimension,
                                                 maxCliplevelInv, faceOffsets, weight);

    if (lowerLevel == upperLevel) return lowSample;

    vec4 highSample = sampleUnifiedClipmapTexture(clipmapTexture, position, upperLevel, voxelSizeL0, volumeDimension,
                                                  maxCliplevelInv, faceOffsets, weight);

    return mix(lowSample, highSample, fract(curLevel));
}

vec4 castConeUnified(in VCTCone c, const VCTScene scene, out VCTIntersection primaryIsect, bool primaryOnly = false) {
    primaryIsect.hit = false;
    // Initialize accumulated color and opacity
    vec4 dst = vec4(0.0);
    // Coefficient used in the computation of the diameter of a cone
    float coneCoefficient = 2.0 * tan(c.aperture * 0.5);

    float startLevel = c.curLevel;
    float voxelSize = scene.voxelSizeL0 * exp2(c.curLevel);

    // Offset startPos in the direction to avoid self occlusion and reduce voxel aliasing
    c.p += c.dir * voxelSize * u_traceDirectionOffset * 0.5;

    float s = 0.0;
    float maxDistance = params.maxDistance;

    // Bound the extent of ray for virtual scene
    if (scene.isVirtual) {
        // 2 is multiplied to avoid missing voxel region for higher clipmap level
        float maxVolumeExtent = scene.voxelSizeL0 * scene.volumeDimension * 2;

        // Find the intersection point with the AABB of virtual objectfor at the cone direction
        // and use the point as a start startLevel for cone casting
        Ray ray = Ray(c.p, c.dir);
        AABBox3D aabb =
            AABBox3D(scene.volumeCenter - maxVolumeExtent * 0.5f, scene.volumeCenter + maxVolumeExtent * 0.5f);
        vec2 tMinMax = rayIntersectsAABB(ray, aabb);
        bool intersected = (tMinMax[0] < tMinMax[1] && tMinMax[0] > 0.0) || inBox(c.p, aabb);
        if (!intersected) return vec4(vec3(0), 1);  // It means no occlusion

        // Start from the object's bounding box if the startPos is not in the box
        // and end at the object's bounding box
        if (tMinMax[0] > 0.f) {
            s = tMinMax[0];
        }
        maxDistance = tMinMax[1];
    }

    // Diameter of cone at start position is the l0 voxel size
    float diameter = max(s * coneCoefficient, scene.voxelSizeL0);

    float stepFactor = max(MIN_STEP_FACTOR, scene.stepFactor);
    float occlusion = 0.0;

    ivec3 faceIndices = computeVoxelFaceIndices(c.dir);  // Implementation in voxelConeTracing/common.glsl
    vec3 faceOffsets = vec3(faceIndices) * FACE_COUNT_INV;
    vec3 weight = c.dir * c.dir;

    float curSegmentLength = voxelSize;

    float minRadius = scene.voxelSizeL0 * scene.volumeDimension * 0.5;

    // Ray marching - compute occlusion and radiance in one go
    while (s < maxDistance && occlusion < 1.0) {
        vec3 position = c.p + c.dir * s;

        float distanceToCenter = length(scene.volumeCenter - position);
        float minLevel = ceil(log2(distanceToCenter / minRadius));

        c.curLevel = log2(diameter / scene.voxelSizeL0);
        c.curLevel = min(max(max(startLevel, c.curLevel), minLevel), scene.maxClipmapLevel - 1);

        // Retrieve radiance by accessing the 3D clipmap (voxel radiance and opacity)
        vec4 radiance =
            sampleUnifiedClipmapLinearly(getRadiance(scene.isVirtual), position, c.curLevel, scene.voxelSizeL0,
                                         scene.volumeDimension, scene.maxClipmapLevelInv, faceOffsets, weight);

        voxelSize = scene.voxelSizeL0 * exp2(c.curLevel);

        // Radiance correction
        float correctionQuotient = curSegmentLength / voxelSize;
        radiance.rgb = radiance.rgb * correctionQuotient;

        // Opacity correction
        float opacity = clamp(1.0 - pow(1.0 - radiance.a, correctionQuotient), 0.0, 1.0);

        vec4 src = vec4(radiance.rgb, opacity);

        // Front-to-back compositing
        dst += clamp(1.0 - dst.a, 0.0, 1.0) * src;
        occlusion += (1.0 - occlusion) * opacity / (1.0 + (s + voxelSize) * params.occlusionDecay);

        // Primary intersection
        if (!primaryIsect.hit && opacity > EPSILON) {
            primaryIsect.position = position;
            primaryIsect.normal = unpackNormal(
                sampleUnifiedClipmapLinearly(getNormal(scene.isVirtual), position, c.curLevel, scene.voxelSizeL0,
                                             scene.volumeDimension, scene.maxClipmapLevelInv, faceOffsets, weight)
                    .rgb);
            primaryIsect.diffuse =
                sampleUnifiedClipmapLinearly(getDiffuse(scene.isVirtual), position, c.curLevel, scene.voxelSizeL0,
                                             scene.volumeDimension, scene.maxClipmapLevelInv, faceOffsets, weight)
                    .rgb;
            primaryIsect.specularA =
                sampleUnifiedClipmapLinearly(getSpecularA(scene.isVirtual), position, c.curLevel, scene.voxelSizeL0,
                                             scene.volumeDimension, scene.maxClipmapLevelInv, faceOffsets, weight);
            primaryIsect.level = c.curLevel;
            primaryIsect.hit = true;
            primaryIsect.t = s;
            if (primaryOnly) break;
        }

        // Step
        float sLast = s;
        s += max(diameter, scene.voxelSizeL0) * stepFactor;
        curSegmentLength = (s - sLast);
        diameter = s * coneCoefficient;
    }

    return clamp(vec4(dst.rgb, 1.0 - occlusion), 0.0, 1.0);
}

// Used to visualize the transition of minLevel selection
vec4 minLevelToColor(float minLevel) {
    vec4 colors[] = {vec4(1.0, 0.0, 0.0, 1.0), vec4(0.0, 1.0, 0.0, 1.0), vec4(0.0, 0.0, 1.0, 1.0),
                     vec4(1.0, 1.0, 0.0, 1.0), vec4(0.0, 1.0, 1.0, 1.0), vec4(1.0, 0.0, 1.0, 1.0),
                     vec4(1.0, 1.0, 1.0, 1.0)};

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

// Diffuse cones for second bounce in real
vec4 castSecondRealDiffuseCones(vec3 startPos, vec3 normal, float realMinLevel) {
    vec4 indirectContribution = vec4(0.0);
    VCTCone c;
    c.p = startPos;  // offset is present inside castConeUnified
    c.curLevel = realMinLevel;
    for (int i = 0; i < DIFFUSE_CONE_COUNT; ++i) {
        VCTIntersection realIsect;
        float cosTheta = dot(normal, DIFFUSE_CONE_DIRECTIONS[i]);

        if (cosTheta < 0.0) continue;

        c.dir = DIFFUSE_CONE_DIRECTIONS[i];
        c.aperture = DIFFUSE_CONE_APERTURE;
        indirectContribution += castConeUnified(c, realScene, realIsect);
    }
    return indirectContribution / (DIFFUSE_CONE_COUNT * 0.5);
}

// Diffuse cones for second bounce in virtual
vec4 castSecondVirtualDiffuseCones(vec3 startPos, vec3 normal, float virtualMinLevel) {
    vec4 indirectContribution = vec4(0.0);
    VCTCone c;
    c.p = startPos;  // offset is present inside castConeUnified
    c.curLevel = virtualMinLevel;
    for (int i = 0; i < DIFFUSE_CONE_COUNT; ++i) {
        VCTIntersection virtualIsect;
        float cosTheta = dot(normal, DIFFUSE_CONE_DIRECTIONS[i]);

        if (cosTheta < 0.0) continue;

        c.dir = DIFFUSE_CONE_DIRECTIONS[i];
        c.aperture = DIFFUSE_CONE_APERTURE;
        indirectContribution += castConeUnified(c, virtualScene, virtualIsect);
    }
    return indirectContribution / (DIFFUSE_CONE_COUNT * 0.5);
}

// Cast diffuse cones at primary intersection point on both real/virtual object
vec4 castDiffuseCones(vec3 startPos, vec3 normal, float realMinLevel, float virtualMinLevel, bool traceVirtual,
                      bool delta = false) {
    vec4 indirectContribution = vec4(0.0);
    VCTCone c;
    c.p = startPos;  // offset is present inside castConeUnified
    for (int i = 0; i < DIFFUSE_CONE_COUNT; ++i) {
        VCTIntersection virtualIsect, realIsect;
        VCTIntersection tmpIsect;
        float cosTheta = dot(normal, DIFFUSE_CONE_DIRECTIONS[i]);

        if (cosTheta < 0.0) continue;

        c.dir = DIFFUSE_CONE_DIRECTIONS[i];
        c.aperture = DIFFUSE_CONE_APERTURE;
        // First hit real
        if (traceVirtual) {
            c.curLevel = virtualMinLevel;
            vec4 virtualRet = castConeUnified(c, virtualScene, virtualIsect) * cosTheta;

            c.curLevel = realMinLevel;
            vec4 realRet = castConeUnified(c, realScene, realIsect) * cosTheta;

            // Second hit virtual
            if (u_secondBounce == 1 && virtualIsect.hit && virtualIsect.t < realIsect.t) {
                // For ambient occlusion (used for DEBUGGING)
                indirectContribution.a += virtualRet.a;

                // Use cone-traced diffuse material of virtual object for smooth result
                vec3 virtualDiffuse = virtualRet.rgb;

                // Virtual object at second bounce (e.g. real->virtual) is assumed to be diffuse
                // So cast diffuse cones
                vec4 secondIntensity =
                    castSecondRealDiffuseCones(virtualIsect.position, virtualIsect.normal, realMinLevel);

                // Radiance from virtual object
                indirectContribution.rgb += secondIntensity.rgb * virtualDiffuse * u_secondIndirectDiffuseIntensity;

                // Anti-radiance from real object (apply occlusion based smoothing)
                indirectContribution.rgb -= (1 - virtualRet.a) * realRet.rgb;
            }
        }
    }
    return indirectContribution / (DIFFUSE_CONE_COUNT * 0.5);
}

void main() {
    float depth = texture2D(u_depthTexture, In.texCoords).r;
    if (depth == 1.0) discard;
    if (u_renderReal == 0 && !isVirtualFrag) discard;
    if (u_renderVirtual == 0 && isVirtualFrag) discard;
    vec3 diffuse = texture(u_diffuseTexture, In.texCoords).rgb;
    vec3 normal = unpackNormal(texture(u_normalMap, In.texCoords).rgb);
    vec3 emission = texture(u_emissionMap, In.texCoords).rgb;
    bool hasEmission = any(greaterThan(emission, vec3(0.0)));
    vec3 posW = worldPosFromDepth(depth);
    vec3 view = normalize(u_eyePos - posW);
    vec4 specColor = texture(u_specularMap, In.texCoords);

    // Scene setup
    virtualScene.voxelSizeL0 = u_virtualVoxelSizeL0;
    virtualScene.volumeDimension = u_virtualVolumeDimension;
    virtualScene.volumeCenter = u_virtualVolumeCenterL0;
    virtualScene.maxClipmapLevel = VIRTUAL_CLIP_LEVEL_COUNT;
    virtualScene.maxClipmapLevelInv = VIRTUAL_CLIP_LEVEL_COUNT_INV;
    virtualScene.stepFactor = u_virtualStepFactor;
    virtualScene.isVirtual = true;

    realScene.voxelSizeL0 = u_voxelSizeL0;
    realScene.volumeDimension = u_volumeDimension;
    realScene.volumeCenter = u_volumeCenterL0;
    realScene.maxClipmapLevel = CLIP_LEVEL_COUNT;
    realScene.maxClipmapLevelInv = CLIP_LEVEL_COUNT_INV;
    realScene.stepFactor = u_stepFactor;
    realScene.isVirtual = false;

    // Parameter setup
    params.occlusionDecay = u_occlusionDecay;
    params.maxDistance = MAX_TRACE_DISTANCE;

    float realMinLevel = getRealMinLevel(posW);
    float virtualMinLevel = getVirtualMinLevel(posW);

    // Compute indirect contribution
    vec4 indirectContribution = vec4(0.0);

    // Offset the starting pos to avoid the self-occlusion
    float voxelSize = frag_voxelSizeL0 * exp2(isVirtualFrag ? virtualMinLevel : realMinLevel);
    vec3 startPosOffset = posW + normal * voxelSize * u_traceStartOffset;
    vec3 startPos = posW;

    if (isVirtualFrag) {
        vec3 virtualIndirectContribution = vec3(0.0);

        float shininess = unpackShininess(specColor.a);
        float roughness = shininessToRoughness(shininess);

        VCTIntersection realIsect, virtualIsect;
        VCTCone c;
        c.aperture = max(roughness, MIN_SPECULAR_APERTURE);
        c.depth = 0;
        const int nSumbsample = 1;

        // Perfect reflection
        if (true) {
            c.dir = reflect(-view, normal);
            c.p = startPos;
            c.aperture = MIN_SPECULAR_APERTURE;

            c.curLevel = realMinLevel;
            vec4 dst = castConeUnified(c, realScene, realIsect);

            c.curLevel = 0;
            c.p = startPosOffset;
            vec4 occlusion = castConeUnified(c, virtualScene, virtualIsect);

            if (virtualIsect.hit && virtualIsect.t < realIsect.t) {
                c.curLevel = virtualIsect.level;
                c.p = virtualIsect.position;
                c.dir = reflect(c.dir, virtualIsect.normal);
                // virtualIndirectContribution.rgb += packNormal(virtualIsect.normal);
                if (u_secondBounce == 1) {
                    // virtualIndirectContribution.rgb += castConeUnified(c, realScene,
                    // realIsect).rgb; virtualIndirectContribution.rgb +=
                    // virtualIsect.diffuse *
                    // castSecondRealDiffuseCones(virtualIsect.position,
                    // virtualIsect.normal, realMinLevel).rgb *
                    // u_virtualIndirectDiffuseIntensity; virtualIndirectContribution.rgb
                    // += packNormal(virtualIsect.normal);
                    virtualIndirectContribution.rgb = occlusion.rgb;
                    // virtualIndirectContribution.rgb = virtualIsect.position;

                    // virtualIndirectContribution.rgb += virtualIsect.diffuse
                    //     * castSecondRealDiffuseCones(virtualIsect.position,
                    //     virtualIsect.normal, realMinLevel).rgb
                    //     * u_virtualIndirectDiffuseIntensity;

                    // virtualIndirectContribution.rgb +=
                    // castSecondVirtualDiffuseCones(virtualIsect.position,
                    // virtualIsect.normal, virtualMinLevel).rgb
                    //     * u_virtualIndirectDiffuseIntensity;
                } else
                    virtualIndirectContribution.rgb += occlusion.rgb;
            } else {
                virtualIndirectContribution.rgb += dst.rgb;
            }
        }
        // Cook-Torrance BRDF
        else
            for (int i = 0; i < nSumbsample; ++i) {
                // sample next direction (for light direction)
                vec2 u = vec2(rand2D(In.texCoords + ivec2(i)), rand2D(In.texCoords + ivec2(i * 2)));
                vec3 direction = sampleBeckmann(u, normal, roughness);
                c.dir = direction;

                c.p = startPos;
                c.curLevel = realMinLevel;
                vec4 lightIntensity = castConeUnified(c, realScene, realIsect);

                c.p = startPosOffset;
                c.curLevel = 0;  // manually set to zero for smoothness
                vec4 occlusion = castConeUnified(c, virtualScene, virtualIsect);

                // 1. Occluded by virtual object
                if (virtualIsect.hit && virtualIsect.t < realIsect.t) {
                    vec4 secondIntensity = vec4(u_ambientSecondIntensity);  // ambient light

                    ivec3 outFaceIndices = computeVoxelFaceIndices(-virtualIsect.normal);
                    vec3 faceOffsets = vec3(outFaceIndices) * FACE_COUNT_INV;
                    vec3 weight = virtualIsect.normal * virtualIsect.normal;
                    vec3 secondDiffuse =
                        sampleUnifiedClipmapLinearly(getDiffuse(true), virtualIsect.position, virtualIsect.level,
                                                     virtualScene.voxelSizeL0, virtualScene.volumeDimension,
                                                     virtualScene.maxClipmapLevelInv, faceOffsets, weight)
                            .rgb;
                    vec4 secondSpecularA = sampleUnifiedClipmapLinearly(
                        getSpecularA(true), virtualIsect.position, virtualIsect.level, virtualScene.voxelSizeL0,
                        virtualScene.volumeDimension, virtualScene.maxClipmapLevelInv, faceOffsets, weight);
                    float shininess2 = unpackShininess(secondSpecularA.a);
                    float roughness2 = shininessToRoughness(shininess2);

                    VCTCone c2;
                    c2.p = virtualIsect.position;
                    c2.curLevel = realMinLevel;
                    // c2.aperture = max(roughness, MIN_SPECULAR_APERTURE);
                    c2.aperture = DIFFUSE_CONE_APERTURE;
                    c2.depth = 0;

                    for (int j = 0; j < 1; ++j) {
                        vec2 u = vec2(rand2D(In.texCoords + ivec2(i + j)), rand2D(In.texCoords + ivec2((i + j) * 2)));
                        vec3 direction = sampleBeckmann(u, virtualIsect.normal, roughness2);
                        c2.dir = direction;
                        vec4 secondIntensity = castConeUnified(c2, realScene, realIsect);

                        virtualIndirectContribution.rgb +=
                            secondIntensity.rgb * diffuse * secondDiffuse * u_virtualIndirectDiffuseIntensity;
                    }
                }
                // 2. Shading using light from real object
                else {
                    vec3 lightDir = -direction;
                    vec3 halfway = normalize(view - lightDir);
                    float nDotL = max(0.0, dot(normal, -lightDir));

                    // Diffuse
                    vec3 F = fresnelSchlick(vec3(0.04), -lightDir, halfway);
                    virtualIndirectContribution += (vec3(1) - F) * nDotL * diffuse * lightIntensity.rgb;

                    // Specular
                    // for the last term F0, use 0.04 if it's plastic, use albeo if it's metal
                    vec3 cook = cookTorranceBRDF(-lightDir, normal, view, halfway, roughness, vec3(0.04));
                    if (any(greaterThan(specColor, vec4(EPSILON))) && any(greaterThan(cook, vec3(EPSILON)))) {
                        // specular
                        virtualIndirectContribution +=
                            cook * specColor.rgb * lightIntensity.rgb * u_indirectSpecularIntensity;
                    }
                }
            }
        virtualIndirectContribution /= nSumbsample;

        indirectContribution.rgb += virtualIndirectContribution;
    } else {
        ivec3 outFaceIndices = computeVoxelFaceIndices(-normal);
        vec3 faceOffsets = vec3(outFaceIndices) * FACE_COUNT_INV;
        vec3 weight = normal * normal;
        // TEMPORAL: Use radiance as a reflectance for now
        vec3 reflectance =
            sampleUnifiedClipmapLinearly(getRadiance(false), startPos, realMinLevel, realScene.voxelSizeL0,
                                         realScene.volumeDimension, realScene.maxClipmapLevelInv, faceOffsets, weight)
                .rgb;

        // Occlusion-based shadow + color bleeding from virtual object
        indirectContribution = castDiffuseCones(startPosOffset, normal, realMinLevel, virtualMinLevel, true, true);

        // if (any(lessThan(indirectContribution.rgb, vec3(0.f))))
        //     indirectContribution.rgb = vec3(1, 0, 0);

        if (u_realReflectance == 1) indirectContribution.rgb *= reflectance;

        indirectContribution.rgb *= u_realIndirectDiffuseIntensity;
    }

    // DIFFUSE_CONE_COUNT includes cones to integrate over a sphere - on the hemisphere there are on average ~half of
    // these cones
    indirectContribution.a *= u_ambientOcclusionFactor;

    // Do not clamp for the anti-radiance
    // indirectContribution = clamp(indirectContribution, 0.0, 1.0);

    vec3 directContribution = vec3(0.0);

    if (hasEmission) {
        directContribution += emission;
    }

    directContribution = clamp(directContribution, 0.0, 1.0);

    out_color = vec4(0.0, 0.0, 0.0, 1.0);

    if ((u_lightingMask & AMBIENT_OCCLUSION_BIT) != 0) directContribution *= indirectContribution.a;

    if ((u_lightingMask & DIRECT_LIGHTING_BIT) != 0) out_color.rgb += directContribution;

    if ((u_lightingMask & INDIRECT_DIFFUSE_LIGHTING_BIT) != 0) out_color.rgb += indirectContribution.rgb;

    // if ((u_lightingMask & INDIRECT_SPECULAR_LIGHTING_BIT) != 0)
    //     out_color.rgb += specularContribution;

    // If only ambient occlusion is selected show ambient occlusion
    if (u_lightingMask == AMBIENT_OCCLUSION_BIT) out_color.rgb = vec3(indirectContribution.a);

    if (u_visualizeMinLevelSelection > 0) {
        if (isVirtualFrag)
            out_color *= minLevelToColor(virtualMinLevel);
        else
            out_color *= minLevelToColor(realMinLevel);
    }

    if (isVirtualFrag) {
        out_color = clamp(out_color, 0.0, 1.0);
    } else {
        // original shading (VCT)
        out_color = clamp(out_color, 0.0, 1.0);
    }
}
