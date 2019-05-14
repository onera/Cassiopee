#version 150 compatibility
layout(triangles) in;
layout(line_strip, max_vertices=2) out;

in Vertex {
    vec4 P0;
    vec4 P1;
    vec4 position;
    vec4 normal;
    vec4 color;
} vertex[];

out vec4 color;
out vec3 Nv;
out vec3 P;
out vec4 vert;
out float gAlpha;

uniform sampler1D colormap;

void main()
{
    int i;
    float ust = (1./3.);
    vec4 bary1 = ust*(vertex[0].P0+vertex[1].P0+vertex[2].P0);
    vec4 bary2 = ust*(vertex[0].P1+vertex[1].P1+vertex[2].P1);
    vec4 bvert= ust*(vertex[0].position+vertex[1].position+vertex[2].position);
    vec4 bcol = ust*(vertex[0].color+vertex[1].color+vertex[2].color);
    float f = length(bcol.rgb);
    f = clamp(f, 0.0f, 1.0f);
    vec3 val = vec3(texture1D(colormap, f));
    bcol = vec4(val.r, val.g, val.b, 1.);
    
    vec4 bnorm = ust*(vertex[0].normal+vertex[1].normal+vertex[2].normal);

    gl_Position = bary1;
    color  = bcol;
    Nv     = bnorm.xyz;
    P      = bary1.xyz;
    vert   = bvert;
    gAlpha = 0.8f;
    EmitVertex();

    gl_Position = bary2;
    color  = bcol;
    Nv     = bnorm.xyz;
    P      = bary1.xyz;
    vert   = bvert;
    gAlpha = 0.f;
    EmitVertex();
    EndPrimitive();
}
