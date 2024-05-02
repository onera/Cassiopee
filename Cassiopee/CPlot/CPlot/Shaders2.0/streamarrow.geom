#version 150 compatibility
layout(triangles) in;
layout(triangle_strip, max_vertices=73) out;

in Vertex {
    vec4 P0;
    vec4 vP;
    vec4 e3;
    vec4 position;
    vec4 normal;
    vec4 color;
} vertex[];

out vec4 color;
out vec3 Nv;
out vec3 P;
uniform int   show_surface;
uniform int   project_vectors;
uniform float density;
uniform sampler1D colormap;
uniform int style_arrow;

void show_trig_surf()
{
    int i;
    // Emission du triangle de base :
    for ( i = 0; i < gl_in.length(); ++i ) {
        gl_Position = vertex[i].P0;
        color  = vec4(0.85,0.85,0.85,1.0);
        Nv     = gl_NormalMatrix * vertex[i].normal.xyz;
        P      = vertex[i].vP.xyz;
        //vert   = vertex[i].position;
        EmitVertex();
    }
    EndPrimitive();
}

float rand(vec2 co){                            //2147483647
    // uniform random
    //return (int(dot(co.xy ,co.xy) * 1664525 + 1013904223)%2147483647)/2147483647.f;
    
    // sinus random
    return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);

    // same with optimal sinus 
    /*
    float x = dot(co.xy ,vec2(12.9898,78.233));
    x *= 0.159155;
    x -= floor(x);
    float xx=x*x;
    float y=-6.87897;
    y=y*xx+33.7755;
    y=y*xx-72.5257;
    y=y*xx+80.5874;
    y=y*xx-41.2408;
    y=y*xx+6.28077;
    return fract(x*y*43758.5453); */
}

void draw_flat_arrow( vec4 origin, vec4 be3, vec4 trn, vec4 bcol )
{
    float ust = (1./3.);
    if ( project_vectors == 1 ) be3 = be3 - dot(be3,trn) * trn;
    vec4 end  = origin + be3;

    vec4 bary1 = gl_ModelViewProjectionMatrix * origin;//ust*(vertex[0].P0+vertex[1].P0+vertex[2].P0);
    vec4 bary2 = gl_ModelViewProjectionMatrix * end;   //ust*(vertex[0].P1+vertex[1].P1+vertex[2].P1);

    trn.xyz = (gl_ModelViewProjectionMatrix * vec4(trn.xyz,0.)).xyz;
    if ( project_vectors == 1 )
    {
        if (trn.z > 0.)
        {
            bary1 -= 1.E-5 * vec4(trn.xyz,0.);
            bary2 -= 1.E-5 * vec4(trn.xyz,0.);
        }
        else
        {
            bary1 += 1.E-5 * vec4(trn.xyz,0.);
            bary2 += 1.E-5 * vec4(trn.xyz,0.);
        }
    }
    be3 = gl_ModelViewProjectionMatrix * be3;
    vec4 be1 = vec4(be3.y, -be3.x, 0., 0.);/*cross(trn.xyz,be3.xyz),0.);*/
    vec3 bnorm = vec3(0., 0., -1.);

    gl_Position = bary1 - 0.125 * be1;
    color  = bcol;
    Nv     = bnorm.xyz;
    P      = bary1.xyz;
    EmitVertex();

    gl_Position = bary1 + 0.125 * be1;
    color  = bcol;
    Nv     = bnorm.xyz;
    P      = bary1.xyz;
    EmitVertex();

    gl_Position = bary1 + 0.75 * be3;//bary2;
    color  = bcol;
    Nv     = bnorm.xyz;
    P      = bary1.xyz;
    EmitVertex();
    EndPrimitive();

    gl_Position = bary1 - 0.2 * be1 + 0.6 * be3;
    color  = bcol;
    Nv     = bnorm.xyz;
    P      = bary1.xyz;
    EmitVertex();

    gl_Position = bary1 + 0.2 * be1 + 0.6 * be3;
    color  = bcol;
    Nv     = bnorm.xyz;
    P      = bary1.xyz;
    EmitVertex();

    gl_Position = bary2;//bary2;
    color  = bcol;
    Nv     = bnorm.xyz;
    P      = bary1.xyz;
    EmitVertex();
    EndPrimitive();
}

void draw_3d_arrow( vec4 origin, vec4 be3, vec4 trn, vec4 bcol )
{
    const float ust = 1./3.;
    // vec4 bnorm= ust*(vertex[0].normal+vertex[1].normal+vertex[2].normal);
    if ( project_vectors == 1 ) be3 = be3 - dot(be3.xyz,trn.xyz) * trn;
    vec4 be1  = vec4(cross(be3.xyz, trn.xyz),0.);
    be3 = gl_ModelViewProjectionMatrix * be3;
    be1 = gl_ModelViewProjectionMatrix * be1;
    float nrmbe1 = length(be1);
    float nrmbe3 = length(be3);
    if ( nrmbe1 < 1.E-1 * nrmbe3 )
    {
        be1 = vec4(be3.y, -be3.x, 0., 0.);
        nrmbe1 = length(be1);
    }
    vec4 vbary = ust*(vertex[0].vP+vertex[1].vP+vertex[2].vP);
    vec4 bary = gl_ModelViewProjectionMatrix * origin;
    trn.xyz = (gl_ModelViewProjectionMatrix * vec4(trn.xyz,0.)).xyz;
    if ( project_vectors == 1 )
    {
        if (trn.z > 0.)
            bary -= 1.E-6 * vec4(trn.xyz,0.);
        else
            bary += 1.E-6 * vec4(trn.xyz,0.);
    }

    gl_Position = bary-0.125*be1;
    color       = bcol;
    Nv          = trn.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();

    gl_Position =  bary+0.125*be1;
    color       = bcol;
    Nv          = trn.xyz;//bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();

    gl_Position = bary + 0.75*be3;
    color       = bcol;
    Nv          = trn.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();
    EndPrimitive();


    gl_Position = bary-0.2*be1 + 0.55 * be3;
    color       = bcol;
    Nv          = trn.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();

    gl_Position =  bary+0.2*be1 + 0.55 * be3;
    color       = bcol;
    Nv          = trn.xyz;//bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();

    gl_Position = bary + be3;
    color       = bcol;
    Nv          = trn.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();
    EndPrimitive();

    // Ortho
    vec4 trn4 = nrmbe1 * normalize(vec4(cross(be1.xyz,be3.xyz),0.));
    /*if ( nrmbe1 > 1.E-6 ) trn4 = normalize(vec4(trn.xyz,0.))*length(be1);
    else trn4 = normalize(vec4(cross(be1.xyz, be3.xyz),0.));*/
    gl_Position = bary-0.125*trn4;
    color       = bcol;
    Nv          = be1.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();

    gl_Position =  bary+0.125*trn4;
    color       = bcol;
    Nv          = be1.xyz;//bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();

    gl_Position = bary + 0.75*be3;
    color       = bcol;
    Nv          = be1.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();
    EndPrimitive();

    gl_Position = bary-0.2*trn4 + 0.55 * be3;
    color       = bcol;
    Nv          = be1.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();

    gl_Position =  bary+0.2*trn4 + 0.55 * be3;
    color       = bcol;
    Nv          = be1.xyz;//bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();

    gl_Position = bary + be3;
    color       = bcol;
    Nv          = be1.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();
    EndPrimitive();    
}

void draw_tetra_arrow( vec4 origin, vec4 be3, vec4 trn, vec4 bcol )
{
    float ust = (1./3.);
    
    if ( project_vectors == 1 ) be3 = be3 - dot(be3,trn) * trn;
    vec4 be1  = vec4(cross(be3.xyz, trn.xyz),0.);
    be3 = gl_ModelViewProjectionMatrix * be3;
    be1 = gl_ModelViewProjectionMatrix * be1;
    float nrmbe1 = length(be1);
    float nrmbe3 = length(be3);
    if ( nrmbe1 < 1.E-1 * nrmbe3 )
    {
        be1 = vec4(be3.y, -be3.x, 0., 0.);
        nrmbe1 = length(be1);
    }
    vec4 vbary= ust*(vertex[0].vP+vertex[1].vP+vertex[2].vP);
    vec4 bary = gl_ModelViewProjectionMatrix * origin;
    trn.xyz = (gl_ModelViewProjectionMatrix * vec4(trn.xyz,0.)).xyz;
    if ( project_vectors == 1 )
    {
        if (trn.z > 0.)
            bary -= 1.E-6 * vec4(trn.xyz,0.);
        else
            bary += 1.E-6 * vec4(trn.xyz,0.);
    }
    vec4 be2 = nrmbe1 * normalize(vec4(cross(be1.xyz,be3.xyz),0.));//length(be3)*vec4(trn.xyz,0.);
    vec3 nbe1 = normalize(be1.xyz).xyz;
    vec3 nbe2 = normalize(be2.xyz).xyz;

    gl_Position = bary-0.15*be1;
    color       = bcol;
    Nv          = nbe2;//bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();

    gl_Position =  bary+0.15*be1;
    color       =  bcol;
    Nv          = nbe2;//bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();

    gl_Position = bary + be3;
    color       = bcol;
    Nv          = nbe2;//bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();
    EndPrimitive();
    // --------------------------------------------------------
    gl_Position = bary-0.15*be2;
    color       = bcol;
    Nv          = nbe1;//bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();

    gl_Position = bary+0.15*be2;
    color       = bcol;
    Nv          = nbe1;//bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();

    gl_Position = bary + be3;
    color       = bcol;
    Nv          = nbe1;//bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();
    EndPrimitive();
}

void generate_uniform_field()
{
    vec4 v0P0 = vertex[0].position;
    vec4 v1P0 = vertex[1].position;
    vec4 v2P0 = vertex[2].position;
    vec4 e1 = v1P0 - v0P0;
    vec4 e2 = v2P0 - v0P0;
    vec4 me1= gl_ModelViewProjectionMatrix * e1;
    vec4 me2= gl_ModelViewProjectionMatrix * e2;

    vec3 unrm = cross(me1.xyz,me2.xyz);

    float nrm = length(unrm);
    float proba = nrm*density;
    int n = int(0.5f*(-1.f+sqrt(1.+8.f*proba)));
    n = ( n > 2 ? 2 : n );
    //int n = int(density);
    if (n>0) {
    for ( int i = 0; i <= n+1; ++i ) {
        float psi = float(i)/(n+2.f);
        for ( int j= 0; j <= n+1-i; ++j ) {
            float ki = float(j)/(n+2.f);
            float te = 1.f-ki-psi;
            vec4 p = psi*v0P0 + ki*v1P0 + te*v2P0;
            vec4 c = psi*vertex[0].color+ki*vertex[1].color+te*vertex[2].color;
            float f = length(c.rgb - vec3(0.5,0.5,0.5))*1.154700543;
            f = clamp(f, 0.0f, 1.0f);
            vec3 val = vec3(texture(colormap, f));
            c = vec4(val.r, val.g, val.b, 1.);

            vec4 nr = psi*vertex[0].normal+ki*vertex[1].normal+te*vertex[2].normal;
            vec4 be3 = psi * vertex[0].e3 + ki * vertex[1].e3 + te * vertex[2].e3;

            if ( style_arrow == 0 ) draw_3d_arrow(p, be3, nr, c);
            else if ( style_arrow == 1 ) draw_flat_arrow(p, be3, nr, c);
            else if ( style_arrow == 2 ) draw_tetra_arrow(p, be3, nr, c);
        }
    } }
    else {
    //if ( n == 0 ) {
        float ust = 0.33333333333;
        // On va tirer une proba pour voir si on emet un vecteur
        vec4 bary1 = ust*(v0P0+v1P0+v2P0);
        vec4 mbary1 = gl_ModelViewProjectionMatrix * bary1;

        //float proba = (nrm_e1+nrm_e2)*0.5f*density;
        float r = 0.5*rand(mbary1.xy);
        if (r < proba)
        {
            vec4 bnorm= ust*(vertex[0].normal+vertex[1].normal+vertex[2].normal);
            vec4 be3 = ust*(vertex[0].e3+vertex[1].e3+vertex[2].e3);
            vec4 bcol = ust*(vertex[0].color+vertex[1].color+vertex[2].color);
            float f = length(bcol.rgb - vec3(0.5,0.5,0.5))*1.154700543;
            f = clamp(f, 0.0f, 1.0f);
            vec3 val = vec3(texture(colormap, f));
            bcol = vec4(val.r, val.g, val.b, 1.);

            if ( style_arrow == 0 ) draw_3d_arrow(bary1, be3, bnorm, bcol);
            else if ( style_arrow == 1 ) draw_flat_arrow(bary1, be3, bnorm, bcol);
            else if ( style_arrow == 2 ) draw_tetra_arrow(bary1, be3, bnorm, bcol);
        }
    }

}

void main()
{
    if ( show_surface == 1 ) show_trig_surf();        // Emission du triangle de base

    // Emission des triangles pour le champs de vecteur
    if ( density < 1.E-6 )
    {
        float ust = (1./3.);
        vec4 bary = ust*(vertex[0].position+vertex[1].position+vertex[2].position);// + vertex[0].translation;
        vec4 be3  = ust*(vertex[0].e3+vertex[1].e3+vertex[2].e3);
        vec4 trn  = normalize(ust*(vertex[0].normal+vertex[1].normal+vertex[2].normal));//cross(be1.xyz,be3.xyz);
        vec4 bcol = ust*(vertex[0].color+vertex[1].color+vertex[2].color);
        float f = length(bcol.rgb - vec3(0.5,0.5,0.5))*1.154700543;
        f = clamp(f, 0.0f, 1.0f);
        vec3 val = vec3(texture(colormap, f));
        bcol = vec4(val.r, val.g, val.b, 1.);
        if ( style_arrow == 0 ) draw_3d_arrow(bary, be3, trn, bcol);
        else if ( style_arrow == 1 ) draw_flat_arrow(bary, be3, trn, bcol);
        else if ( style_arrow == 2 ) draw_tetra_arrow(bary, be3, trn, bcol);
    }
    else generate_uniform_field();
}
