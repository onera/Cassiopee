#version 400 compatibility
//
// Anisotropic light (brushed metal)
//
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
/*varying vec3 MCposition;
varying vec4 color;
varying vec3 normal;
varying vec3 tangent;
varying vec3 light;
varying vec3 eye;
varying vec4 vertex;*/
uniform float Scale;

void main(void) 
{
    //MCposition = vec3(gl_Vertex) * Scale;
    v2f_out.vdata1 = gl_Vertex * Scale;
    vec3 t = vec3(1.,0.1,.1);
    t = normalize(t);
    //normal = normalize(gl_NormalMatrix * gl_Normal);
    v2f_out.view_normal = vec4(gl_NormalMatrix * gl_Normal,0.);
    v2f_out.nrm_view_normal = normalize(v2f_out.view_normal);

    //vec3 t = vec3(normal.z,normal.z,-normal.x-normal.y);
    //t = normalize(t);

    // tangent = t - dot(t,normal)*normal; // strand dir
    // tangent = normalize(tangent);
    v2f_out.vdata2 = normalize(vec4(t - dot(t,v2f_out.view_normal.xyz)*v2f_out.view_normal.xyz, 0.));

    //if (dot(tangent, tangent) < 1.e-12)
    //{
    //  t = vec3(0,1.,0.);
    //  tangent = t - dot(t,normal)*normal;
    //}   
    //bitangent = normalize(cross(tangent, normal));

    //vec3 P = vec3(gl_ModelViewMatrix * gl_Vertex);
    //eye = -P;
    v2f_out.mv_position = gl_ModelViewMatrix * gl_Vertex;
    // vec3 lightPosition = gl_LightSource[0].position.xyz;
    // light = lightPosition - P;
    v2f_out.vdata3      = gl_LightSource[0].position - v2f_out.mv_position;
    // vertex = gl_Vertex;
    v2f_out.position = gl_Vertex;
    // color = gl_Color;
    v2f_out.color = gl_Color;

    v2ct_out.position = gl_Vertex;
    v2ct_out.color    = gl_Color ;
    v2ct_out.data_comp= ivec4(1,1,1,0);
    v2ct_out.vdata1 = v2f_out.vdata1;
    v2ct_out.vdata2 = v2f_out.vdata2;
    v2ct_out.vdata3 = v2f_out.vdata3;

    gl_Position = ftransform();
}
