// Procedural bricks
varying vec2 MCposition;
varying vec3 Nv;
varying vec3 P;
varying vec4 initColor;
varying vec4 vertex;

void main(void)
{
    P = vec3(gl_ModelViewMatrix * gl_Vertex);
    Nv = gl_NormalMatrix * gl_Normal;
    vertex = gl_Vertex;
    MCposition = gl_Vertex.xy;
    initColor = gl_Color;
    gl_Position = ftransform();
}