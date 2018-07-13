#version 150 compatibility

out Vertex {
    vec4 P0;
    vec4 P1;
    vec4 position;
    vec4 normal;
    vec4 color;
} vertex;
uniform int fix_length;
uniform float scale;

void main()
{
    gl_Position = ftransform();
    vertex.position = gl_Vertex;
    vertex.P0 = gl_Position;
    vec4 transpos = vec4(gl_Color.xyz,0.) - vec4(0.5,0.5,0.5,0.0);
    transpos = (1-fix_length)*transpos + fix_length*normalize(transpos);

    vertex.P1 = gl_ModelViewProjectionMatrix*(gl_Vertex+scale*vec4(transpos.xyz,0.));
    vertex.normal   = vec4(gl_Normal,0.);
    vertex.color    = gl_Color;
}
