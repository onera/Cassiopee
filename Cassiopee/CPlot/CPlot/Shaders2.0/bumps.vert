//
// Vertex shader for procedural bumps
//
varying vec3 LightDir;
varying vec3 EyeDir;
varying vec4 color;

void main() 
{
    EyeDir         = vec3(gl_ModelViewMatrix * gl_Vertex);
    gl_Position    = ftransform();
    gl_TexCoord[0] = gl_MultiTexCoord0;
    
    vec3 n = normalize(gl_NormalMatrix * gl_Normal);
    vec3 Tangent = vec3(n.z,n.z,-n.x-n.y);
    Tangent = normalize(Tangent);
    vec3 t = normalize(gl_NormalMatrix * Tangent);
    vec3 b = cross(n, t);

    vec3 v;
    vec3 LightPosition = gl_LightSource[0].position.xyz;
    v.x = dot(LightPosition, t);
    v.y = dot(LightPosition, b);
    v.z = dot(LightPosition, n);
    LightDir = normalize(v);

    v.x = dot(EyeDir, t);
    v.y = dot(EyeDir, b);
    v.z = dot(EyeDir, n);
    EyeDir = normalize(v);
    color = gl_Color;
}