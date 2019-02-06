#version 400 compatibility
/* XRay shader */ 

in V2F_OUT
{
    vec4 position;
    vec4 mv_position;
    vec4 mvp_position;
    vec4 view_normal;
    vec4 nrm_view_normal;
    vec4 color;
    vec4 vdata1, vdata2, vdata3, vdata4;
} v2f_out;

/*varying vec3 N;
varying vec3 I;
varying vec4 color;*/

uniform float EdgeFalloff;
uniform float intensity;

void main()
{
    vec3 N = v2f_out.view_normal.xyz;
    vec3 I = v2f_out.mv_position.xyz;
    vec4 color = v2f_out.color;
    float opacity = dot(normalize(N), normalize(I));
    opacity = abs(opacity);
    opacity = 1.0 - pow(opacity, EdgeFalloff);
    opacity = clamp(opacity, 0., 1.);
    
    //vec4 result = pow(color, vec4(1.0 / 2.2));
    gl_FragColor = color*(1.+opacity*2.*(intensity-1.));
    gl_FragColor.a = opacity;
}
