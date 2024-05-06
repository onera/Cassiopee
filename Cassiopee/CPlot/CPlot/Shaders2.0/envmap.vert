// envmap shader (chrome)
varying vec3  Normal;
varying vec3  EyeDir;
varying vec4 color;
varying vec4 vertex;

void main(void) 
{
    Normal = normalize(gl_NormalMatrix * gl_Normal);
    vec4 pos = gl_ModelViewMatrix * gl_Vertex;
    EyeDir = pos.xyz;
    vertex = gl_Vertex;
    color = gl_Color;
    gl_Position = ftransform();
}