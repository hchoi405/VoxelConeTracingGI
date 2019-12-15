#version 430
#extension GL_ARB_shading_language_include : enable

#include "/voxelConeTracing/voxelizationGeom.glsl"

layout(points) in;
layout(points, max_vertices = 1) out;

in Vertex
{
    vec3 normalW;
    vec2 uv;
    vec3 color;
} In[1];

out Geometry
{
	vec3 posW;
    vec3 normalW;
    vec2 uv;
    vec3 color;
} Out;

void main()
{
    // int idx = getDominantAxisIdx(gl_in[0].gl_Position.xyz, gl_in[1].gl_Position.xyz, gl_in[2].gl_Position.xyz);
    int idx = 0;
	gl_ViewportIndex = idx; 
	
	for (int i = 0; i < 1; ++i)
    {
        gl_Position = u_viewProj[idx] * gl_in[i].gl_Position;
        Out.uv = In[i].uv;
		Out.posW = gl_in[i].gl_Position.xyz;
        Out.normalW = In[i].normalW;
        Out.color = In[i].color;
        
        EmitVertex();
    }

    EndPrimitive();
}
