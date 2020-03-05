#ifndef COMMON_STRUCT_H
#define COMMON_STRUCT_H

struct VCTIntersection {
    bool hit;
    vec3 position;
    vec3 normal;
    float level;
    float t;
};

struct VCTCone {
    vec3 p;
    vec3 dir;
    float aperture;
    float curLevel;
    int depth;
};

struct VCTScene {
    float voxelSizeL0;
    uint volumeDimension;
    vec3 volumeCenter;
    uint maxClipmapLevel;
    float maxClipmapLevelInv;
    bool isVirtual;
    float stepFactor;
};

struct VCTParams {
    float occlusionDecay;
    float maxDistance;
};

#endif