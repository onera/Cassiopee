#version 150 compatibility
layout(triangles) in;
layout(triangle_strip, max_vertices=15) out;

in Vertex {
    vec4 P0;
    vec4 vP;
    vec4 e1;
    vec4 e3;
    vec4 position;
    vec4 normal;
    vec4 color;
    vec4 translation;
} vertex[];

out vec4 color;
out vec3 Nv;
out vec3 P;
out vec4 vert;
out float gAlpha;
uniform int   show_surface;

uniform sampler1D colormap;

void main()
{
    int i;
    if ( show_surface == 1 ) {
        // Emission du triangle de base :
        for ( i = 0; i < gl_in.length(); ++i ) {
	       gl_Position = vertex[i].P0;
	       color  = vec4(vertex[i].color.xyz,0.5); //vec4(0.5,0.35,0.5,1.0);
           Nv     = vertex[i].normal.xyz;
	       P      = vertex[i].vP.xyz;
	       vert   = vertex[i].position;
	       gAlpha = 0.f;
	       EmitVertex();
        }
        EndPrimitive();
    }
    // Emission des triangles pour le champs de vecteur
    float ust = (1./3.);
    vec4 bary = ust*(vertex[0].P0+vertex[1].P0+vertex[2].P0) + vertex[0].translation;
    vec4 vbary = ust*(vertex[0].vP+vertex[1].vP+vertex[2].vP);
    vec4 be1  = ust*(vertex[0].e1+vertex[1].e1+vertex[2].e1);
    vec4 be3  = ust*(vertex[0].e3+vertex[1].e3+vertex[2].e3);
    vec4 bcol = vec4(ust*(vertex[0].color+vertex[1].color+vertex[2].color).xyz,1.);
    float f = length(bcol.rgb - vec3(0.5,0.5,0.5))*2./1.73205080;
    f = clamp(f, 0.0f, 1.0f);
    vec3 val = vec3(texture(colormap, f));
    bcol = vec4(val.r, val.g, val.b, 1.);
    
    vec4 bvert= ust*(vertex[0].position+vertex[1].position+vertex[2].position);
    // vec4 bnorm= ust*(vertex[0].normal+vertex[1].normal+vertex[2].normal);
    vec3 trn  = cross(be1.xyz,be3.xyz);

    gl_Position = bary-0.125*be1;
    color       = bcol;
    Nv          = trn.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    vert        = bvert;
    gAlpha      = 0.f;
    EmitVertex();

    gl_Position =  bary+0.125*be1;
    color       = bcol;
    Nv          = trn.xyz;//bnorm.xyz;
    P           = vbary.xyz;
    vert        = bvert;
    gAlpha      = 0.f;
    EmitVertex();

    gl_Position = bary + 0.75*be3;
    color       = bcol;
    Nv          = trn.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    vert        = bvert;
    gAlpha      = 0.f;
    EmitVertex();
    EndPrimitive();


    gl_Position = bary-0.3*be1 + 0.55 * be3;
    color       = bcol;
    Nv          = trn.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    vert        = bvert;
    gAlpha      = 0.0f;
    EmitVertex();

    gl_Position =  bary+0.3*be1 + 0.55 * be3;
    color       = bcol;
    Nv          = trn.xyz;//bnorm.xyz;
    P           = vbary.xyz;
    vert        = bvert;
    gAlpha      = 0.0f;
    EmitVertex();

    gl_Position = bary + be3;
    color       = bcol;
    Nv          = trn.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    vert        = bvert;
    gAlpha      = 0.f;
    EmitVertex();
    EndPrimitive();


    // Ortho
    vec4 trn4 = normalize(vec4(trn,0.))*length(be1);
    gl_Position = bary-0.125*trn4;
    color       = bcol;
    Nv          = be1.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    vert        = bvert;
    gAlpha      = 0.f;
    EmitVertex();

    gl_Position =  bary+0.125*trn4;
    color       = bcol;
    Nv          = be1.xyz;//bnorm.xyz;
    P           = vbary.xyz;
    vert        = bvert;
    gAlpha      = 0.f;
    EmitVertex();

    gl_Position = bary + 0.75*be3;
    color       = bcol;
    Nv          = be1.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    vert        = bvert;
    gAlpha      = 0.f;
    EmitVertex();
    EndPrimitive();

    gl_Position = bary-0.3*trn4 + 0.55 * be3;
    color       = bcol;
    Nv          = be1.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    vert        = bvert;
    gAlpha      = 0.0f;
    EmitVertex();

    gl_Position =  bary+0.3*trn4 + 0.55 * be3;
    color       = bcol;
    Nv          = be1.xyz;//bnorm.xyz;
    P           = vbary.xyz;
    vert        = bvert;
    gAlpha      = 0.0f;
    EmitVertex();

    gl_Position = bary + be3;
    color       = bcol;
    Nv          = be1.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    vert        = bvert;
    gAlpha = 0.f;
    EmitVertex();
    EndPrimitive();
}
