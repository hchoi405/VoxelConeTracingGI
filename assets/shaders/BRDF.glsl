#ifndef BRDF_GLSL
#define BRDF_GLSL

#define BLINN_PHONG_MODE_IDX 0
#define COOK_TORRANCE_MODE_IDX 1

#ifndef PI
#define PI 3.14159265358979323846
#define PI_2 1.57079632679489661923   // PI/2
#define PI_4 0.785398163397448309616  // PI/4
#endif

#ifndef INV_PI
#define INV_PI 0.318309886183790671538  // 1 / PI
#define INV_2PI 0.15915494309189533577  // 1 / 2*PI
#endif

#include "/util/noise.glsl"

// etaI: commonly the vacuum (1.0)
// etaT: the glass (1.5-1.6)
float fresnelDieletric(vec3 view, vec3 normal, float etaI, float etaT)
{
    float cosThetaI = dot(view, normal);

    // Normal should be directed to view direction
    if (cosThetaI < 0) return 0;
    
    float sinThetaI = sqrt(max(0, 1 - cosThetaI * cosThetaI));
    float sinThetaT = etaI / etaT * sinThetaI;

    // Total internal reflection
    if (sinThetaT >= 1) return 1;
    float cosThetaT = sqrt(max(0, 1 - sinThetaT * sinThetaT));
    float Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
                  ((etaT * cosThetaI) + (etaI * cosThetaT));
    float Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
                  ((etaI * cosThetaI) + (etaT * cosThetaT));
    return (Rparl * Rparl + Rperp * Rperp) / 2;
}

// Calculate the refracted ray over 
// etaRatio: etaI / etaT
bool refract(const vec3 view, const vec3 n, float etaRatio, out vec3 wt) {
    // Compute $\cos \theta_\roman{t}$ using Snell's law
    float cosThetaI = dot(n, view);
    float sin2ThetaI = max(0, 1 - cosThetaI * cosThetaI);
    float sin2ThetaT = etaRatio * etaRatio * sin2ThetaI;

    // Handle total internal reflection for transmission
    if (sin2ThetaT >= 1) return false;
    float cosThetaT = sqrt(1 - sin2ThetaT);
    wt = etaRatio * -view + (etaRatio * cosThetaI - cosThetaT) * n;
    return true;
}

float fresnel(vec3 viewVector, vec3 halfway)
{
    return pow(1.0 - dot(viewVector, halfway), 5);
}

vec3 fresnelSchlick(vec3 F0, vec3 l, vec3 n)
{
    float p = pow(1.0 - max(0.0, dot(l, n)), 5);
    return F0 + (1.0 - F0) * p;
}

vec3 fresnelSchlickMicrofacet(vec3 F0, vec3 l, vec3 h)
{
    float p = pow(1.0 - dot(l, h), 5);
    return F0 + (1.0 - F0) * p;
}

vec3 blinnPhongBRDF(vec3 Id, vec3 kd, vec3 Is, vec3 ks, vec3 normal, vec3 lightVec, vec3 halfway, float shininess)
{
    float diffuseFactor = max(dot(normal, lightVec), 0.0);
    vec3 diffuse = Id * kd * diffuseFactor;
    
    float specFactor = pow(max(dot(normal, halfway), 0.0), shininess);
    vec3 specular = vec3(0.0);
    
    if (shininess > 0.0 && diffuseFactor > 0.0)
        specular = Is * ks * specFactor;
    
    return diffuse + specular;
}

vec3 sampleBeckmann(vec2 u, vec3 n, float roughness)
{
    float phi = u.x * 2 * PI;
    float theta = acos(1.0 / (1 - roughness*roughness*log(1-u.y)));
    float sinTheta = sin(theta);
    float cosPhi = cos(phi);
    float sinPhi = sin(phi);
    float cosTheta = cos(theta);
    vec3 ph = vec3(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);

    // Make arbitrary coordinate using normal
    vec3 tangent, bitangent;
    if (abs(n.x) > abs(n.y))
        tangent = vec3(-n.z, 0, n.x) / sqrt(n.x * n.x + n.z * n.z);
    else
        tangent = vec3(0, n.z, -n.y) / sqrt(n.y * n.y + n.z * n.z);
    bitangent = cross(n, tangent);

    /* Make our hemisphere orient around the normal. */
    return tangent * ph.x + bitangent * ph.y + n * ph.z;
}

// return: sampled point at local coordinate
vec2 concentricSampleDisk(vec2 u) {
    // Map uniform random numbers to $[-1,1]^2$
    vec2 uOffset = 2.f * u - vec2(1, 1);

    // Handle degeneracy at the origin
    if (uOffset.x == 0 && uOffset.y == 0) return vec2(0, 0);

    // Apply concentric mapping to point
    float theta, r;
    if (abs(uOffset.x) > abs(uOffset.y)) {
        r = uOffset.x;
        theta = PI_4 * (uOffset.y / uOffset.x);
    } else {
        r = uOffset.y;
        theta = PI_2 - PI_4 * (uOffset.x / uOffset.y);
    }
    return r * vec2(cos(theta), sin(theta));
}

// return: sampled direction at local coordinate
vec3 sampleCosineHemisphere(vec2 u) {
    vec2 d = concentricSampleDisk(u);
    const float z = sqrt(max(0.f, 1 - dot(d, d)));
    return vec3(d.x, d.y, z);
}

vec3 sampleUniformHemisphere(vec2 u, vec3 n) {
    float z = u[0];
    float r = sqrt(max(0.f, 1.f - z * z));
    float phi = 2 * PI * u[1];
    float theta = acos(z);
    // return Vector3f(r * std::cos(phi), r * std::sin(phi), z);

    float sinTheta = sin(theta);
    float cosPhi = cos(phi);
    float sinPhi = sin(phi);
    float cosTheta = cos(theta);
    vec3 ph = vec3(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);

    // Make arbitrary coordinate using normal
    vec3 tangent, bitangent;
    if (abs(n.x) > abs(n.y))
        tangent = vec3(-n.z, 0, n.x) / sqrt(n.x * n.x + n.z * n.z);
    else
        tangent = vec3(0, n.z, -n.y) / sqrt(n.y * n.y + n.z * n.z);
    bitangent = cross(n, tangent);
    
    /* Make our hemisphere orient around the normal. */
    return tangent * ph.x + bitangent * ph.y + n * ph.z;

}

float beckmannNDF(vec3 n, vec3 h, float roughness)
{
    float rSq = roughness * roughness;
    float nh = max(0.0, dot(n,h));
    float nhSq = nh * nh;
    float nhPow4 = nhSq * nhSq;
    
    float c0 = 1.0 / (rSq * nhPow4);
    float c1 = exp(-(1.0 - nhSq) / (nhSq * rSq));
    
    return c0 * c1;
}

vec3 cookTorranceBRDF(vec3 l, vec3 n, vec3 v, vec3 h, float roughness, vec3 F0)
{
    vec3 F = fresnelSchlick(F0, l, h);
    
    float nh = max(0.0, dot(n,h));
    float vh = max(0.0, dot(v,h));
    float nv = max(0.0, dot(n,v));
    float nl = max(0.0, dot(n,l));
    float G0 = (2.0 * nh * nv) / vh;
    float G1 = (2.0 * nh * nl) / vh;
    
    float G = min(1.0, max(0.0, min(G0, G1)));
    
    float D = beckmannNDF(n, h, roughness);
    
    vec3 numerator = F * G * D;
    float denominator = 4.0 * nl * nv;
    
    return max(vec3(0.0), numerator / denominator);
}

float mean(vec3 v) { return (v.x + v.y + v.z) / 3.f; }

bool isSameHemisphere(const vec3 w1, const vec3 w2) {
  return w1.z * w2.z > 0;
}

float absCosTheta(const vec3 w) { return abs(w.z); }

vec3 reflect2(const vec3 w) { return vec3(-w.x, -w.y, w.z); }

// Returns the probability of the transfer of energy between two direction
// (local)
float phongPdf(const vec3 w_out, const vec3 w_in, float kappa, float exponent) {
    if (!isSameHemisphere(w_in, w_out)) return 0;

    // FIXME Diffusive lobe sampling
    const float pdf_d = absCosTheta(w_out) * INV_PI;

    // Specular lobe sampling
    const float alpha = dot(w_out, reflect2(w_in));
    const float pdf_s = alpha > 0 ? pow(alpha, exponent) * (exponent + 1.0f) / (2.0f * PI) : 0;

    return kappa * pdf_s + (1 - kappa) * pdf_d;
}

// Return sampled incident direction
vec3 phongSample(const vec3 w_out, uint seed, float kappa, float exponent, out bool sampleSpecular) {
    vec2 u = rand2D(seed);

    sampleSpecular = false;
    if (u.x < kappa) {
        u.x /= kappa;
    } else {
        u.x = (u.x - kappa) / (1 - kappa);
        sampleSpecular = true;
    }

    if (sampleSpecular) {
        float sinAlpha = sqrt(1.0f - pow(u.y, 2.0f / (exponent + 1)));
        float cosAlpha = pow(u.y, 1.0f / (exponent + 1));
        float phi = 2.0f * PI * u.x;
        return vec3(sinAlpha * cos(phi), sinAlpha * sin(phi), cosAlpha);
    }

    return sampleCosineHemisphere(u);
}

// Evaluate brdf
vec3 phongEval(const vec3 w_out, const vec3 w_in, vec3 Kd, vec3 Ks, float exponent) {
    if (!isSameHemisphere(w_in, w_out)) return vec3(0.f);

    vec3 result = vec3(0.f);

    // Diffusive part
    result += Kd * INV_PI;

    // Specular part
    float alpha = dot(w_out, reflect2(w_in));
    if (alpha > 0) result += Ks * ((exponent + 2) * INV_2PI * pow(alpha, exponent));

    return result * absCosTheta(w_out);
}

// PBRT-v3 like sample_f
// return in local coordinate
vec3 phongSample_f(vec3 w_out, out vec3 w_in, uint seed, out float p, 
                    vec3 Kd, vec3 Ks, float exponent, out bool sampleSpecular) {    
    vec2 u = rand2D(seed);
    const float avg_diffuse = mean(Kd);
    const float avg_specular = mean(Ks);
    const float kappa = avg_specular / (avg_diffuse + avg_specular);

    w_in = phongSample(w_out, seed, kappa, exponent, sampleSpecular);
    if (w_out.z < 0) w_in.z *= -1;
    p = phongPdf(w_out, w_in, kappa, exponent);
    if (p == 0) return vec3(0.f);
    return phongEval(w_out, w_in, Kd, Ks, exponent);
}

/**
 * Generate a cosine-weighted random point on the unit hemisphere oriented around 'n'.
 *
 * @param rand a vector containing pseudo-random values
 * @param n    the normal to orient the hemisphere around
 * @returns    the cosine-weighted point on the oriented hemisphere
 */
vec3 randomCosineWeightedHemispherePoint(vec3 rand, vec3 n, out float pdf) {
    float r = rand.x * 0.5 + 0.5;       // [-1..1) -> [0..1)
    float angle = (rand.y + 1.0) * PI;  // [-1..1] -> [0..2*PI)
    float sr = sqrt(r);
    vec2 p = vec2(sr * cos(angle), sr * sin(angle));
    /*
     * Unproject disk point up onto hemisphere:
     * 1.0 == sqrt(x*x + y*y + z*z) -> z = sqrt(1.0 - x*x - y*y)
     */
    vec3 ph = vec3(p.xy, sqrt(1.0 - p * p));

    // PDF
    pdf = ph.z * INV_PI;

    /*
     * Compute some arbitrary tangent space for orienting
     * our hemisphere 'ph' around the normal. We use the camera's up vector
     * to have some fix reference vector over the whole screen.
     */
    vec3 tangent = normalize(rand);
    vec3 bitangent = cross(tangent, n);
    tangent = cross(bitangent, n);

    /* Make our hemisphere orient around the normal. */
    return tangent * ph.x + bitangent * ph.y + n * ph.z;
}

#endif