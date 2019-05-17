#version 400 compatibility
//
// Vertex shader for producing clouds (mostly sunny)
//
uniform float Scale;

out V2F_OUT
{
    vec4 position;
    vec4 mv_position;
    vec4 mvp_position;
    vec4 view_normal;
    vec4 nrm_view_normal;
    vec4 color;
    vec4 vdata1, vdata2, vdata3, vdata4;
} v2f_out;

out V2CT_OUT
{
    vec4 position;
    vec4 color;
    ivec4 data_comp; // Raconte si vdata1,2,3 ou 4 est utilise
    vec4 vdata1, vdata2, vdata3, vdata4;
} v2ct_out;

/*varying float LightIntensity;
varying vec3  MCposition;
varying vec4 initColor;
varying vec4 vertex;*/

void main(void)
{
    vec3 LightPos = gl_LightSource[0].position.xyz;
    vec3 ECposition = vec3(gl_ModelViewMatrix * gl_Vertex);
    vec3 L = LightPos - ECposition;
    vec3 Nv = gl_NormalMatrix * gl_Normal;
    vec3 N = normalize(Nv);
    if (dot(N, L) < 0.) N = -N;
    float LightIntensity  = dot(normalize(L), N)*1.5;
    //MCposition = vec3(gl_Vertex) * Scale;
    v2f_out.vdata1 = vec4(vec3(gl_Vertex) * Scale, LightIntensity);
    // initColor
    v2f_out.color = gl_Color;
    //vertex = gl_Vertex;
    v2f_out.position = gl_Vertex;

    v2ct_out.position = gl_Vertex;
    v2ct_out.color    = gl_Color ;
    v2ct_out.data_comp= ivec4(1,0,0,0);
    v2ct_out.vdata1   = vec4(vec3(gl_Vertex) * Scale, LightIntensity);
    
    gl_Position = ftransform();
}
