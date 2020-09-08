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
const float MIN_STEP_FACTOR = 0.01;
// const float MIN_SPECULAR_APERTURE = 0.1; // 5.73 degrees
float MIN_SPECULAR_APERTURE = 0.05;  // 2.86 degrees
// const float MIN_SPECULAR_APERTURE = 0.025;  // 2.86 degrees

uniform sampler2D u_diffuseTexture;
uniform sampler2D u_normalMap;
uniform sampler2D u_specularMap;
uniform sampler2D u_emissionMap;
uniform sampler2D u_depthTexture;
uniform sampler2D u_virtualMap;
uniform sampler2D u_backgroundTexture;

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
uniform vec3 u_virtualMin;
uniform vec3 u_virtualMax;

uniform int u_BRDFMode;
uniform mat4 u_viewProjInv;
uniform uint u_volumeDimension;
uniform vec3 u_eyePos;
uniform float u_voxelSizeL0;
uniform vec3 u_volumeCenterL0;
uniform float u_stepFactor;
uniform float u_virtualStepFactor;
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
uniform uint u_toggleViewCone;
uniform float u_viewAperture;
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
uniform int u_extraStep;
uniform float u_glassEta;
uniform float u_phongShininess;
uniform int u_rotateCone;
uniform float u_localRatio;
uniform int u_excludeEmptyFace;
uniform int u_material;
uniform int u_subsample;

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
#else  

// 16 cones for lower quality (8 on average per hemisphere)
const vec3 DIFFUSE_CONE_NORMAL = vec3(0, 1, 0);

const int DIFFUSE_CONE_COUNT = 16;
const float DIFFUSE_CONE_APERTURE = 0.931517;  // 50 degree
const vec3 DIFFUSE_CONE_DIRECTIONS[16] = {
    vec3(0.57735, 0.57735, 0.57735),       vec3(0.57735, -0.57735, -0.57735),     vec3(-0.57735, 0.57735, -0.57735),
    vec3(-0.57735, -0.57735, 0.57735),     vec3(-0.903007, -0.182696, -0.388844), vec3(-0.903007, 0.182696, 0.388844),
    vec3(0.903007, -0.182696, 0.388844),   vec3(0.903007, 0.182696, -0.388844),   vec3(-0.388844, -0.903007, -0.182696),
    vec3(0.388844, -0.903007, 0.182696),   vec3(0.388844, 0.903007, -0.182696),   vec3(-0.388844, 0.903007, 0.182696),
    vec3(-0.182696, -0.388844, -0.903007), vec3(0.182696, 0.388844, -0.903007),   vec3(-0.182696, 0.388844, 0.903007),
    vec3(0.182696, -0.388844, 0.903007)};

// // For hemisphere sampling
// const int DIFFUSE_CONE_COUNT = 8;
// const float DIFFUSE_CONE_APERTURE = 0.931517;
// const vec3 DIFFUSE_CONE_DIRECTIONS[8] = {
//     vec3(0.57735, 0.57735, 0.57735),
//     vec3(-0.57735, 0.57735, -0.57735),
//     vec3(-0.903007, 0.182696, 0.388844),
//     vec3(0.903007, 0.182696, -0.388844),
//     vec3(0.388844, 0.903007, -0.182696),
//     vec3(-0.388844, 0.903007, 0.182696),
//     vec3(0.182696, 0.388844, -0.903007),
//     vec3(-0.182696, 0.388844, 0.903007)
// };

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

bool inBox(vec3 p, vec3 bmin, vec3 bmax) { return all(greaterThan(p, bmin)) && all(lessThan(p, bmax)); }

bool inBox(vec3 p, AABBox3D box) { return all(greaterThan(p, box.minPos)) && all(lessThan(p, box.maxPos)); }

bool inSphere(vec3 p, vec3 center, float radius) { return length(p - center) < radius; }

bool inLocal(vec3 p, vec3 min, vec3 max, float ratio) {
    vec3 extent = max - min;
    return inBox(p, min - ratio * extent, max + ratio * extent);
}

float dx[8] = {0, 0, 0, 0, 1, 1, 1, 1};
float dy[8] = {0, 0, 1, 1, 0, 0, 1, 1};
float dz[8] = {0, 1, 0, 1, 0, 1, 0, 1};

vec4 interpolateTrilinearly(bool neighborExist[8], vec4 neighborColor[8], vec3 frac, out bool exist) {
    bool bv01, bv23, bv0123, bv45, bv67, bv4567, bv01234567;
    vec4 v01, v23, v0123, v45, v67, v4567, v01234567;

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

    bv01234567 = true;
    if (bv0123 && bv4567) {
        v01234567 = mix(v0123, v4567, frac.x);
    } else if (bv0123) {
        v01234567 = v0123;
    } else if (bv4567) {
        v01234567 = v4567;
    } else {
        v01234567 = vec4(0.0);
    }

    exist = bv01234567;
    if (bv01234567)
        return v01234567;
    else
        return vec4(0);
}

// All input texture should have GL_NEAREST filter to properly work
vec4 sampleClipmapTexture2(sampler3D clipmapTexture, vec3 posW, int clipmapLevel, float voxelSizeL0,
                           uint volumeDimension, float maxCliplevelInv, vec3 faceOffsets, vec3 weight) {
    float voxelSize = voxelSizeL0 * exp2(clipmapLevel);
    float extent = voxelSize * volumeDimension;

    vec3 neighborPos[8];
    bool neighborExist[3][8];
    vec4 neighborColor[3][8];

    vec3 frac = fract(posW / voxelSize);

    for (uint i = 0; i < 8; ++i) {
        neighborPos[i] = posW + vec3(dx[i], dy[i], dz[i]) * voxelSize;

#ifdef VOXEL_TEXTURE_WITH_BORDER
        neighborPos[i] = (fract(neighborPos[i] / extent) * volumeDimension + vec3(BORDER_WIDTH)) /
                         (float(volumeDimension) + 2.0 * BORDER_WIDTH);
#else
        neighborPos[i] = fract(neighborPos[i] / extent);
#endif
        neighborPos[i].y += clipmapLevel;
        neighborPos[i].y *= maxCliplevelInv;
        neighborPos[i].x *= FACE_COUNT_INV;

        for (uint j = 0; j < 3; ++j) {
            neighborExist[j][i] = false;
            vec4 rgba = texture(clipmapTexture, neighborPos[i] + vec3(faceOffsets[j], 0.0, 0.0));
            if (rgba.a > EPSILON) {
                neighborExist[j][i] = true;
                neighborColor[j][i] = rgba;
            }
        }
    }

    vec4 sum = vec4(0);
    vec4 dirSum[3];
    bool dirExist[3] = {false, false, false};
    float weightSum = 0.0;
    for (int i = 0; i < 3; ++i) {
        dirSum[i] = interpolateTrilinearly(neighborExist[i], neighborColor[i], frac, dirExist[i]);
        if (dirExist[i]) {
            weightSum += weight[i];
            sum += dirSum[i] * weight[i];
            // return clamp(vec4(vec3(weight[i]), 1), 0.0, 1.0);
        }
    }
    if (weightSum > 0.0) sum /= weightSum;

    return clamp(sum, 0.0, 1.0);
}

vec4 sampleClipmapLinearly2(sampler3D clipmapTexture, vec3 position, float curLevel, float voxelSizeL0,
                            uint volumeDimension, float maxCliplevelInv, vec3 faceOffsets, vec3 weight) {
    int lowerLevel = int(floor(curLevel));
    int upperLevel = int(ceil(curLevel));

    vec4 lowSample = sampleClipmapTexture2(clipmapTexture, position, lowerLevel, voxelSizeL0, volumeDimension,
                                           maxCliplevelInv, faceOffsets, weight);

    if (lowerLevel == upperLevel) return lowSample;

    vec4 highSample = sampleClipmapTexture2(clipmapTexture, position, upperLevel, voxelSizeL0, volumeDimension,
                                            maxCliplevelInv, faceOffsets, weight);

    return mix(lowSample, highSample, fract(curLevel));
}

vec4 sampleClipmapTexture(sampler3D virtualClipmapTexture, vec3 posW, int clipmapLevel, float voxelSizeL0,
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

vec4 sampleClipmapLinearly(sampler3D clipmapTexture, vec3 position, float curLevel, float voxelSizeL0,
                           uint volumeDimension, float maxCliplevelInv, vec3 faceOffsets, vec3 weight) {
    int lowerLevel = int(floor(curLevel));
    int upperLevel = int(ceil(curLevel));

    vec4 lowSample = sampleClipmapTexture(clipmapTexture, position, lowerLevel, voxelSizeL0, volumeDimension,
                                          maxCliplevelInv, faceOffsets, weight);

    if (lowerLevel == upperLevel) return lowSample;

    vec4 highSample = sampleClipmapTexture(clipmapTexture, position, upperLevel, voxelSizeL0, volumeDimension,
                                           maxCliplevelInv, faceOffsets, weight);

    return mix(lowSample, highSample, fract(curLevel));
}

vec4 castCone(in VCTCone c, const VCTScene scene, out VCTIntersection isect) {
    // Initialize the intersection
    isect.hit = false;
    isect.normal = vec3(0.0);
    isect.diffuse = vec3(0.0);
    isect.specularA = vec4(0.0);

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

    bool startCount = false;
    int terminateCounter = u_extraStep;

    // Ray marching - compute occlusion and radiance in one go
    while (s < maxDistance && occlusion < 1.0) {
        vec3 position = c.p + c.dir * s;

        float distanceToCenter = length(scene.volumeCenter - position);
        float minLevel = ceil(log2(distanceToCenter / minRadius));

        c.curLevel = log2(diameter / scene.voxelSizeL0);
        c.curLevel = min(max(max(startLevel, c.curLevel), minLevel), scene.maxClipmapLevel - 1);
        voxelSize = scene.voxelSizeL0 * exp2(c.curLevel);

        // Retrieve radiance by accessing the 3D clipmap (voxel radiance and opacity)
        vec4 radiance = sampleClipmapLinearly(getRadiance(scene.isVirtual), position, c.curLevel, scene.voxelSizeL0,
                                              scene.volumeDimension, scene.maxClipmapLevelInv, faceOffsets, weight);

        // Radiance correction (diffuse, specular)
        float correctionQuotient = curSegmentLength / voxelSize;
        radiance.rgb = radiance.rgb * correctionQuotient;

        // Opacity correction
        float opacity = clamp(1.0 - pow(1.0 - radiance.a, correctionQuotient), 0.0, 1.0);

        // Sample only opaque voxel
        if (radiance.a > EPSILON) {
            vec3 normal =
                unpackNormal(sampleClipmapLinearly(getNormal(scene.isVirtual), position, c.curLevel, scene.voxelSizeL0,
                                                   scene.volumeDimension, scene.maxClipmapLevelInv, faceOffsets, weight)
                                 .rgb);
            vec3 diffuse = sampleClipmapLinearly(getDiffuse(scene.isVirtual), position, c.curLevel, scene.voxelSizeL0,
                                                 scene.volumeDimension, scene.maxClipmapLevelInv, faceOffsets, weight)
                               .rgb;
            vec4 specularA =
                sampleClipmapLinearly(getSpecularA(scene.isVirtual), position, c.curLevel, scene.voxelSizeL0,
                                      scene.volumeDimension, scene.maxClipmapLevelInv, faceOffsets, weight);

            normal = normal * correctionQuotient;
            diffuse = diffuse * correctionQuotient;
            specularA = specularA * correctionQuotient;

            isect.normal += clamp(1.0 - dst.a, 0.0, 1.0) * normal;
            isect.diffuse += clamp(1.0 - dst.a, 0.0, 1.0) * diffuse;
            isect.specularA += clamp(1.0 - dst.a, 0.0, 1.0) * specularA;
        }

        // Front-to-back compositing
        // Primary intersection
        if (!isect.hit && opacity > EPSILON) {
            isect.position = position;
            isect.level = c.curLevel;
            isect.hit = true;
            isect.t = s;
            startCount = true;
        }
        vec4 src = vec4(radiance.rgb, opacity);
        dst += clamp(1.0 - dst.a, 0.0, 1.0) * src;
        occlusion += (1.0 - occlusion) * opacity / (1.0 + (s + voxelSize) * params.occlusionDecay);

        // Step
        float sLast = s;
        s += max(diameter, scene.voxelSizeL0) * stepFactor;
        curSegmentLength = (s - sLast);
        diameter = s * coneCoefficient;

        if (startCount) {
            // if (terminateCounter-- < 0) break;
        }
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
vec4 castSecondDiffuseConesToReal(vec3 startPos, vec3 normal, float realMinLevel) {
    vec4 indirectContribution = vec4(0.0);
    VCTCone c;
    c.p = startPos;  // offset is present inside castCone
    c.curLevel = realMinLevel;
    for (int i = 0; i < DIFFUSE_CONE_COUNT; ++i) {
        VCTIntersection realIsect;
        float cosTheta = dot(normal, DIFFUSE_CONE_DIRECTIONS[i]);

        if (cosTheta < 0.0) continue;

        c.dir = DIFFUSE_CONE_DIRECTIONS[i];
        c.aperture = DIFFUSE_CONE_APERTURE;
        indirectContribution += castCone(c, realScene, realIsect) * cosTheta;
    }
    return indirectContribution / (DIFFUSE_CONE_COUNT * 0.5);
}

// Diffuse cones for second bounce in virtual
vec4 castSecondDiffuseConesToVirtual(vec3 startPos, vec3 normal, float virtualMinLevel) {
    vec4 indirectContribution = vec4(0.0);
    VCTCone c;
    c.p = startPos;  // offset is present inside castCone
    c.curLevel = virtualMinLevel;
    for (int i = 0; i < DIFFUSE_CONE_COUNT; ++i) {
        VCTIntersection virtualIsect;
        float cosTheta = dot(normal, DIFFUSE_CONE_DIRECTIONS[i]);

        if (cosTheta < 0.0) continue;

        c.dir = DIFFUSE_CONE_DIRECTIONS[i];
        c.aperture = DIFFUSE_CONE_APERTURE;
        indirectContribution += castCone(c, virtualScene, virtualIsect);
    }
    return indirectContribution / (DIFFUSE_CONE_COUNT * 0.5);
}

mat3 skewMatrix(vec3 v) {
    return mat3(0, v[2], -v[1], -v[2], 0, v[0], v[1], -v[0], 0);
}

// Cast diffuse cones at primary intersection point on both real/virtual object
vec4 castDiffuseCones(vec3 startPos, vec3 normal, float realMinLevel, float virtualMinLevel, bool traceVirtual, uint seed = 0, bool rotateCone = false) {
    vec4 indirectContribution = vec4(0.0);
    VCTCone c;
    c.p = startPos;  // offset is present inside castCone

    // Random rotaiton on sphere (from Fast Random Rotation Matrices [Arvo 91])
    mat3 M = mat3(1.f);

    if (rotateCone) {
        float x1 = rand(seed), x2 = rand(seed), x3 = rand(seed);
        float sqrt3 = sqrt(x3);
        vec3 v = vec3(cos(PI2 * x2) * sqrt3, sin(PI2 * x2), sqrt(1 - x3));
        float PI2x1 = PI2 * x1;

        // Column major (first 3 elements for first column)
        mat3 R = mat3(
            cos(PI2x1),-sin(PI2x1), 0,
            sin(PI2x1), cos(PI2x1), 0,
            0, 0, 1
        );
        mat3 v2 = mat3(
            v[0] * v[0], v[0] * v[1], v[0] * v[2],
            v[1] * v[0], v[1] * v[1], v[1] * v[2],
            v[2] * v[0], v[2] * v[1], v[2] * v[2]
        );
        mat3 H = mat3(1.f) - 2 * v2;
        M = -H * R;
    }

    // // Alignment of cones to normal
    // mat3 R = mat3(1.f);
    // if (u_debugFlag == 1) {
    //     vec3 v = cross(DIFFUSE_CONE_NORMAL, normal);
    //     float s = length(v);
    //     float cc = dot(DIFFUSE_CONE_NORMAL, normal);
    //     mat3 vs = skewMatrix(v);
    //     R = mat3(1.f) + vs + vs * vs * (1 - cc) / (s * s);
    // }

    for (int i = 0; i < DIFFUSE_CONE_COUNT; ++i) {
        VCTIntersection virtualIsect, realIsect;
        VCTIntersection tmpIsect;
        vec3 newDir = M * DIFFUSE_CONE_DIRECTIONS[i];
        // vec3 newDir = R * DIFFUSE_CONE_DIRECTIONS[i];
        float cosTheta = dot(normal, newDir);

        if (cosTheta < 0.0f) continue;
        // if (cosTheta < 0.0f) {
        //     return vec4(1, 0, 0, 1);
        // }

        c.dir = newDir;
        c.aperture = DIFFUSE_CONE_APERTURE;
        // First hit real
        if (traceVirtual) {
            // c.curLevel = virtualMinLevel;
            c.curLevel = 0;  // Use zero level to remove sphere artifacts
            vec4 virtualRet = castCone(c, virtualScene, virtualIsect) * cosTheta;

            c.curLevel = realMinLevel;
            vec4 realRet = castCone(c, realScene, realIsect) * cosTheta;

            // Second hit virtual
            if (virtualIsect.hit && virtualIsect.t < realIsect.t) {
                // For ambient occlusion (used for DEBUGGING)
                indirectContribution.a += virtualRet.a;

                // Use cone-traced diffuse material of virtual object for smooth result
                vec3 virtualDiffuse = virtualRet.rgb;

                // Virtual object at second bounce (e.g. real->virtual) is assumed to be diffuse
                // So cast diffuse cones
                vec4 secondIntensity = vec4(1);
                if (u_secondBounce == 1)
                    secondIntensity =
                        castSecondDiffuseConesToReal(virtualIsect.position, normalize(virtualIsect.normal), realMinLevel);

                // Radiance from virtual object
                indirectContribution.rgb += secondIntensity.rgb * virtualDiffuse * u_secondIndirectDiffuseIntensity;

                // Anti-radiance from real object (apply occlusion based smoothing)
                indirectContribution.rgb -= (1 - virtualRet.a) * realRet.rgb;
            }
        }
        // First hit virtual
        else {
            c.curLevel = realMinLevel;
            vec4 realRet = castCone(c, realScene, realIsect) * cosTheta;

            c.curLevel = 0;
            vec4 virtualRet = castCone(c, virtualScene, virtualIsect) * cosTheta;

            if (u_secondBounce == 1 && virtualIsect.hit && virtualIsect.t < realIsect.t) {
                vec3 secondNormal = normalize(virtualIsect.normal);
                vec3 secondDiffuse = virtualIsect.diffuse;  // intersected color (artifacts)
                //     vec3(0.880392f, 0.768627f, 0.323725f);  // For self-intersection in rise103 diffuse Buddha
                // vec3(1, 0, 0);  // For inter-intersection in dasan106 two buddha

                vec4 secondIntensity =
                    castSecondDiffuseConesToReal(virtualIsect.position, secondNormal, realMinLevel) * cosTheta;

                // Radiance from virtual object (second bounce)
                indirectContribution.rgb +=
                    // Virtual-virtual occlusion
                    (1 - virtualRet.a)
                    // Second bounce contribution (assume diffuse)
                    * secondDiffuse * secondIntensity.rgb * u_secondIndirectDiffuseIntensity;

                // Radiance from real object (first bounce)
                indirectContribution.rgb += virtualRet.a * realRet.rgb;
            } else {
                indirectContribution.rgb += realRet.rgb * cosTheta;
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
    const vec3 normal = unpackNormal(texture(u_normalMap, In.texCoords).rgb);
    vec3 emission = texture(u_emissionMap, In.texCoords).rgb;
    vec3 background = texture(u_backgroundTexture, In.texCoords).rgb;
    bool hasEmission = any(greaterThan(emission, vec3(0.0)));
    vec3 posW = worldPosFromDepth(depth);
    vec3 view = normalize(posW - u_eyePos);
    vec4 specColor = texture(u_specularMap, In.texCoords);

    // Check whether it's local
    bool isLocal = inLocal(posW, u_virtualMin, u_virtualMax, u_localRatio);
    MIN_SPECULAR_APERTURE = u_viewAperture;

    ShadingFrame frame;
    frame.n = normal;
    coordinateSystem(frame.n, frame.b, frame.t);
    vec3 w_out = frame.to_local(-view);

    // DIFFERENTIAL RENDERING (mask)
    // if (isVirtualFrag) out_color = vec4(1);
    // else out_color = vec4(0);
    // return;

    // DEBUG: Sample at posW
    // ivec3 faceIndices = computeVoxelFaceIndices(-normal);
    // vec3 faceOffsets = vec3(faceIndices) * FACE_COUNT_INV;
    // // Real weight
    // vec3 w = normal * normal;
    // out_color = sampleClipmapLinearly2(getRadiance(false), posW, getRealMinLevel(posW), u_voxelSizeL0,
    //                                    u_volumeDimension, CLIP_LEVEL_COUNT_INV, faceOffsets, w);
    // return;

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

    vec3 directContribution = vec3(0.0);

    // Non-local region
    if (!isLocal) {
        out_color = vec4(background, 1);
        return;
    }

    // Parameter setup
    params.occlusionDecay = u_occlusionDecay;
    params.maxDistance = MAX_TRACE_DISTANCE;

    float realMinLevel = getRealMinLevel(posW);
    float virtualMinLevel = getVirtualMinLevel(posW);

    // DEBUG: Visualize elvel
    // if (isVirtualFrag) {
    //     out_color = vec4(vec3(virtualMinLevel), 1);
    //     return;
    // }
    // else {
    //     out_color = vec4(vec3(realMinLevel), 1);
    //     return;
    // }

    // DEBUG: Cast cone from u_eyePos
    if (u_toggleViewCone == 1) {
        VCTIntersection tmpIsect;
        VCTCone tmpCone;
        tmpCone.depth = 0;
        tmpCone.dir = -view;
        tmpCone.p = u_eyePos;
        tmpCone.aperture = u_viewAperture;
        tmpCone.curLevel = getRealMinLevel(u_eyePos);
        out_color.rgb = castCone(tmpCone, realScene, tmpIsect).rgb;
        return;
    }

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

        // Perfect reflection (Mirror)
        if (u_material == 0) {
        // if (diffuse[0] > 0.49f && diffuse[0] < 0.51f) { // Mirror ball
            c.dir = reflect(view, normal);
            c.p = startPos;
            c.aperture = MIN_SPECULAR_APERTURE;
            c.curLevel = realMinLevel;
            vec4 dst = castCone(c, realScene, realIsect);

            c.curLevel = virtualMinLevel;
            c.p = startPosOffset;
            vec4 occlusion = castCone(c, virtualScene, virtualIsect);

            // Virtual-Virtual occlusion
            if (u_secondBounce == 1 && virtualIsect.hit && virtualIsect.t < realIsect.t) {
                c.curLevel = virtualIsect.level;
                c.p =
                    virtualIsect.position + normalize(virtualIsect.normal) * u_virtualVoxelSizeL0 * u_traceStartOffset;
                c.dir = reflect(c.dir, normalize(virtualIsect.normal));

                vec4 real = castCone(c, realScene, realIsect);
                vec4 virtual2 = castCone(c, virtualScene, virtualIsect);

                if (virtualIsect.hit && virtualIsect.t < realIsect.t) {
                    vec3 secondNormal = normalize(virtualIsect.normal);
                    vec3 secondDiffuse = virtualIsect.diffuse;
                    // vec3 secondDiffuse = vec3(0.48f, 0.392f, 0.114f); // Yellow glossy lucy dasan613
                    vec4 secondIntensity = 
                        castSecondDiffuseConesToReal(virtualIsect.position, secondNormal, realMinLevel);

                    // Radiance from virtual object
                    virtualIndirectContribution.rgb +=
                        secondIntensity.rgb * secondDiffuse * u_secondIndirectDiffuseIntensity;
                } else {
                    virtualIndirectContribution.rgb += real.rgb;
                }
            } else {
                virtualIndirectContribution.rgb += dst.rgb;
            }
            virtualIndirectContribution.rgb *= u_indirectSpecularIntensity;
        }
        // Glass
        else if (u_material == 1) {
            float etaI = 1.0;                                       // Incident medium (Vacuum)
            float etaT = u_glassEta;                                // Transmitted medium (Glass, 1.5-1.6)
            vec3 refractedDir;
            if (refract(view, normal, etaI / etaT, refractedDir)) {
                // Go further to prevent self occlusion
                c.p = startPos - normal * voxelSize * u_traceStartOffset;
                c.aperture = MIN_SPECULAR_APERTURE;
                c.curLevel = virtualMinLevel;
                c.dir = refractedDir;

                castCone(c, virtualScene, virtualIsect);
                // Inside glass object
                if (virtualIsect.hit) {
                    // // DEBUG: visualize the normal of inside
                    // virtualIndirectContribution = packNormal(normalize(virtualIsect.normal));

                    // Switched eta
                    if (refract(c.dir, normalize(virtualIsect.normal), etaI / etaT, refractedDir)) {
                        // Outside glass object
                        c.curLevel = realMinLevel;
                        c.p = virtualIsect.position;
                        c.curLevel = getRealMinLevel(u_eyePos);  // Use this due to some artifacts
                        c.dir = refractedDir;
                        virtualIndirectContribution = castCone(c, realScene, realIsect).rgb;
                        // virtualIndirectContribution = vec3(0, 1, 0);
                    } else {
                        // total reflrection (false color)
                        virtualIndirectContribution = vec3(0, 0, 1);
                    }
                }
                // Wrong hit check or thin glass
                else {
                    virtualIndirectContribution = vec3(1, 0, 0);
                }
            } else {
                // total reflrection (false color)
                virtualIndirectContribution = vec3(0, 0, 1);
            }
            virtualIndirectContribution *= u_indirectSpecularIntensity;
        }
        // Diffuse
        else if (u_material == 2) {
        // else if (diffuse[0] > 0.8f) { // Red Bunny
            for (int i = 0; i < u_subsample; ++i) {
                uint seed = seed(uint(gl_FragCoord.y) * 640 + uint(gl_FragCoord.x), i);
                virtualIndirectContribution +=
                    castDiffuseCones(startPosOffset, normal, realMinLevel, virtualMinLevel, false, seed).rgb;
                if (u_rotateCone == 0) break; // No need for subsample without rotating the cone
            }
            if (u_rotateCone == 1) virtualIndirectContribution /= u_subsample;
            virtualIndirectContribution *= diffuse;

            virtualIndirectContribution *= u_virtualIndirectDiffuseIntensity;
        }
        // Phong
        else if (u_material == 3) {
        // else if (diffuse[0] < 0.499f) { // Yellow lucy
            for (int i = 0; i < u_subsample; ++i) {
                // sample next direction (for light direction)
                uint see = seed(uint(gl_FragCoord.y) * 640 + uint(gl_FragCoord.x), i);
                vec3 Kd = diffuse;
                vec3 Ks = vec3(1.f);
                vec3 sumK = Kd + Ks;
                float maxFactor = max(max(sumK[0], sumK[1]), sumK[2]);
                if (maxFactor > 1.f) {
                    Kd = 0.99f * Kd / maxFactor;
                    Ks = 0.99f * Ks / maxFactor;
                }
                float exponent = u_phongShininess;
                float pdf;
                vec3 direction;

                bool sampleSpecular = false;
                vec3 phongBRDF = phongSample_f(w_out, direction, see, pdf, Kd, Ks, exponent, sampleSpecular);

                direction = frame.to_world(direction);
                c.dir = direction;
                c.p = startPosOffset;

                if (sampleSpecular) {
                    vec3 specularContribution = vec3(0);
                    c.aperture = MIN_SPECULAR_APERTURE;

                    c.curLevel = realMinLevel;
                    vec4 realRet = castCone(c, realScene, realIsect);

                    c.curLevel = 0;  // manually set to zero for smoothness
                    vec4 virtualRet = castCone(c, virtualScene, virtualIsect);

                    // Occluded by virtual object
                    if (u_secondBounce == 1 && virtualIsect.hit && virtualIsect.t < realIsect.t) {
                        vec3 secondNormal = normalize(virtualIsect.normal);
                        vec3 secondDiffuse = virtualIsect.diffuse;  // intersected color (artifacts)

                        vec4 secondIntensity =
                            castSecondDiffuseConesToReal(virtualIsect.position, secondNormal, realMinLevel);

                        // Radiance from virtual object (second bounce)
                        specularContribution.rgb +=
                            // Virtual-virtual occlusion
                            (1 - virtualRet.a)
                            // Second bounce contribution (assume diffuse)
                            * secondDiffuse * secondIntensity.rgb * u_secondIndirectDiffuseIntensity;

                        // Radiance from real object (first bounce)
                        specularContribution.rgb += virtualRet.a * realRet.rgb;
                    }
                    // Shading using light from real object
                    else {
                        specularContribution += realRet.rgb * phongBRDF / pdf;
                    }
                    virtualIndirectContribution += specularContribution * u_indirectSpecularIntensity;
                } else {
                    virtualIndirectContribution +=
                        u_virtualIndirectDiffuseIntensity *
                        castDiffuseCones(startPosOffset, normal, realMinLevel, virtualMinLevel, false, see).rgb *
                        phongBRDF / pdf;
                }
            }
            virtualIndirectContribution /= u_subsample;
        }

        indirectContribution.rgb += virtualIndirectContribution;
    } else {
        // out_color.rgb = packNormal(normal);
        // return;
        ivec3 outFaceIndices = computeVoxelFaceIndices(-normal);
        vec3 faceOffsets = vec3(outFaceIndices) * FACE_COUNT_INV;
        vec3 weight = normal * normal;
        vec3 reflectance = vec3(0.263517f, 0.23897f, 0.22877);  // From reflectance estimation (neon)

        for (int i = 0; i < u_subsample; ++i) {
            uint seed = seed(uint(gl_FragCoord.y) * 640 + uint(gl_FragCoord.x), i);
            // Occlusion-based shadow + color bleeding from virtual object
            indirectContribution += castDiffuseCones(startPosOffset, normal, realMinLevel, virtualMinLevel, true, seed, true);
        }
        indirectContribution /= u_subsample;

        if (u_realReflectance == 1) indirectContribution.rgb *= reflectance;

        indirectContribution.rgb *= u_realIndirectDiffuseIntensity;

        indirectContribution.rgb += background;
    }

    // DIFFUSE_CONE_COUNT includes cones to integrate over a sphere - on the hemisphere there are on average ~half of
    // these cones
    indirectContribution.a *= u_ambientOcclusionFactor;

    // Do not clamp for the anti-radiance
    // indirectContribution = clamp(indirectContribution, 0.0, 1.0);

    if (hasEmission) {
        directContribution += emission;
    }

    directContribution = clamp(directContribution, 0.0, 1.0);

    out_color = vec4(0.0, 0.0, 0.0, 1.0);

    // if ((u_lightingMask & AMBIENT_OCCLUSION_BIT) != 0) directContribution *= (1.f - indirectContribution.a);

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

    out_color = clamp(out_color, 0.0, 1.0);
}
