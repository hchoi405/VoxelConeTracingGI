#ifndef BRDF_GLSL
#define BRDF_GLSL

#define BLINN_PHONG_MODE_IDX 0
#define COOK_TORRANCE_MODE_IDX 1

#ifndef PI
#define PI 3.14159265
#endif

#ifndef INV_PI
#define INV_PI 0.318309886183791
#endif

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

/**
 * Generate a cosine-weighted random point on the unit hemisphere oriented around 'n'.
 * 
 * @param rand a vector containing pseudo-random values
 * @param n    the normal to orient the hemisphere around
 * @returns    the cosine-weighted point on the oriented hemisphere
 */
vec3 randomCosineWeightedHemispherePoint(vec3 rand, vec3 n, out float pdf) {
    float r = rand.x * 0.5 + 0.5;      // [-1..1) -> [0..1)
    float angle = (rand.y + 1.0) * PI; // [-1..1] -> [0..2*PI)
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