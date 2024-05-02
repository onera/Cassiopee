#version 400 compatibility
// Gooch (hand drawing) shader
/*varying float NdotL;
varying vec3 ReflectVec;
varying vec3 ViewVec;
varying vec3 ecPos;
varying vec3 tnorm;
varying vec4 color;
varying vec4 vertex;*/
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

void main(void)
{
    vec3 ecPos = vec3(gl_ModelViewMatrix * gl_Vertex);
    v2f_out.mv_position = vec4(ecPos,1.);   
    vec3 tnorm = normalize(gl_NormalMatrix * gl_Normal);
    v2f_out.view_normal = vec4(tnorm,0.);
    v2f_out.nrm_view_normal = vec4(normalize(tnorm),0);
    vec3 lightVec = normalize(gl_LightSource[0].position.xyz - ecPos);
    //ReflectVec = normalize(reflect(-lightVec, tnorm));
    v2f_out.vdata1 = vec4(normalize(reflect(-lightVec, tnorm)),0);
    //ViewVec = normalize(-ecPos);
    float NdotL = (dot(lightVec, tnorm) + 1.0) * 0.5;
    v2f_out.vdata2 = vec4(normalize(-ecPos), NdotL);
    /*vertex = gl_Vertex;
    color = gl_Color;*/
    v2f_out.position = gl_Vertex;
    v2f_out.color    = gl_Color;

    v2ct_out.position = gl_Vertex;
    v2ct_out.color    = gl_Color ;
    v2ct_out.data_comp= ivec4(1,1,0,0);
    v2ct_out.vdata1   = vec4(normalize(reflect(-lightVec, tnorm)),0);
    v2ct_out.vdata2   = vec4(normalize(-ecPos.xyz), NdotL);

    gl_Position = ftransform();
}
