#version 150 compatibility

out Vertex {
    vec4 P0;
    vec4 P1;
    vec4 normal;
    vec4 translation;
    vec4 color;
} vertex;

uniform int fix_length;
uniform float scale;

void main()
{
    gl_Position = ftransform();
    vertex.P0    = gl_Vertex;
    vertex.color = gl_Color;
    vec4 transpos = vec4(gl_Color.xyz,0.) - vec4(0.5,0.5,0.5,0.0);
    transpos = (1-fix_length)*transpos + 0.5*fix_length*normalize(transpos);
    vertex.P1 = vertex.P0 + scale * transpos;
    vertex.translation = gl_ModelViewProjectionMatrix * scale * transpos;
}
