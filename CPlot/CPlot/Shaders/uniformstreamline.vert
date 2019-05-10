#version 150 compatibility

out Vertex {
    vec4 P0;
    vec4 translation;
    vec4 position;
    vec4 normal;
    vec4 color;
} vertex;
uniform float scale;
uniform int fix_length;

void main()
{
    gl_Position = ftransform();
    vertex.position = gl_Vertex;
    vertex.P0 = gl_Position;
    vec4 transpos = vec4(gl_Color.xyz,0) - vec4(0.5,0.5,0.5,0.0);
    transpos = (1-fix_length)*transpos + fix_length*normalize(transpos);
    vertex.translation = gl_ModelViewProjectionMatrix*(scale*vec4(transpos.xyz,0.));
    //vertex.translation = (1-fix_length)*vertex.translation + fix_length*normalize(vertex.translation);
    vertex.normal = vec4(gl_Normal,0.);
    vertex.color  = gl_Color;
}
