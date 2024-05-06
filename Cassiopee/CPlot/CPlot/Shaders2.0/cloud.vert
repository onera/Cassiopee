//
// Vertex shader for producing clouds (mostly sunny)
//
uniform float Scale;
varying float LightIntensity;
varying vec3  MCposition;
varying vec4 initColor;
varying vec4 vertex;

void main(void)
{
    vec3 LightPos = gl_LightSource[0].position.xyz;
    vec3 ECposition = vec3(gl_ModelViewMatrix * gl_Vertex);
    vec3 L = LightPos - ECposition;
    MCposition = vec3(gl_Vertex) * Scale;
    vec3 Nv = gl_NormalMatrix * gl_Normal;
    vec3 N = normalize(Nv);
    if (dot(N, L) < 0.) N = -N;
    LightIntensity  = dot(normalize(L), N)*1.5;
    initColor = gl_Color;
    vertex = gl_Vertex;
    gl_Position = ftransform();
}
