#ifndef NOISE_GLSL
#define NOISE_GLSL

uint tea(uint val0, uint val1) {
    uint v0 = val0;
    uint v1 = val1;
    uint s0 = 0;

    for (uint n = 0; n < 4; n++) {
        s0 += 0x9e3779b9;
        v0 += ((v1 << 4) + 0xa341316c) ^ (v1 + s0) ^ ((v1 >> 5) + 0xc8013ea4);
        v1 += ((v0 << 4) + 0xad90777d) ^ (v0 + s0) ^ ((v0 >> 5) + 0x7e95761e);
    }

    return v0;
}

uint seed(uint pixel_index, uint frame_index) { return tea(pixel_index, frame_index); }

// Generate random uint in [0, 2^24)
uint lcg(uint prev) {
    const uint LCG_A = 1664525u;
    const uint LCG_C = 1013904223u;
    prev = (LCG_A * prev + LCG_C);
    return prev & 0x00FFFFFF;
}

// Generate random float in [0, 1)
float rand(uint prev) { return float(lcg(prev)) / float(0x01000000); }

// Generate two random float in [0, 1)
vec2 rand2D(uint prev) { return vec2(rand(prev), rand(prev * prev)); }
vec2 rand2D(uvec2 prev) { return vec2(rand(prev.x), rand(prev.y)); }

#endif