#version 410 compatibility

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

// Boolean indiquant si la surface doit etre representee ou non.
uniform int show_surface;
// Boolean racontant si les vecteurs doivent etre projetes ou non sur la surface.
uniform int project_vectors;
// Grille cartesienne servant a generer une densite de vecteurs donnee.
uniform float density;
// La palette de couleur utilisee pour representer le champs
uniform sampler1D colormap;
// Le style de fleche utilise pour representer les vecteurs
uniform int style_arrow;

/* But : afficher la surface portant le champs de vecteur
*/
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
int seed;

/*
But:
----
    Retourner un reel unitaire pseudo aleatoire avec une loi uniforme

Discussion:
-----------
    Cette routine met en oeuvre la recursion suivante :

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    L'arithmetique entiere ne demande jamais plus que 32 bits, bit de signe compris

Reference:
----------
    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

Parametres:
-----------
    Input/output, int *seed, une graine pour le generateur de nombre aleatoire

    Output, un reel simple precision avec la valeur pseudo-aleatoire
*/
float uniform_rand()
{
    const int i4_huge = 2147483647;
    int k = seed / 127773;
    seed = 16807 * ( seed - k * 127773 ) - k * 2836;
    if ( seed < 0 ) seed = seed + i4_huge;
    return abs(fract(sin(seed*12.9898/2147483647.) * 43758.5453));
/*    const int i4_huge = 2147483647;
    int k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 ) seed = seed + i4_huge;
    return float(( seed ) * 4.656612875E-10);*/
}

float rand(vec3 co){                            //2147483647
    // uniform random
    //return (int(dot(co.xy ,co.xy) * 1664525 + 1013904223)%2147483647)/2147483647.f;
    
    // sinus random
    return fract(sin(dot(co.xyz ,vec3(12.9898,78.233, 83.453))) * 43758.5453);

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
    vec4 dir = normalize(be3);
    vec4 be1 = vec4(dir.y, -dir.x, 0., 0.);/*cross(trn.xyz,be3.xyz),0.);*/
    vec3 bnorm = vec3(0., 0., -1.);
    const float thickness = 0.01f;
    const float head_arrow_width = 1.5;
    vec4 offset = vec4(be1 * bary1.w * thickness / 2.0f );

    gl_Position = bary1 - offset;
    color  = bcol;
    Nv     = bnorm.xyz;
    P      = bary1.xyz;
    EmitVertex();

    gl_Position = bary1 + offset;
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

    gl_Position = bary1 - head_arrow_width * offset + 0.6 * be3;
    color  = bcol;
    Nv     = bnorm.xyz;
    P      = bary1.xyz;
    EmitVertex();

    gl_Position = bary1 + head_arrow_width * offset + 0.6 * be3;
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
    be1 = normalize(be1);
    vec4 dir = normalize(be3);
    const float thickness = 0.01f;
    const float head_arrow_width = 1.5;
    vec4 offset = vec4(be1 * bary.w * thickness / 2.0f );

    gl_Position = bary-offset;
    color       = bcol;
    Nv          = trn.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();

    gl_Position =  bary+offset;
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


    gl_Position = bary- head_arrow_width * offset + 0.55 * be3;
    color       = bcol;
    Nv          = trn.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();

    gl_Position =  bary+ head_arrow_width * offset + 0.55 * be3;
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
    vec4 trn4 = /*nrmbe1 **/ normalize(vec4(cross(be1.xyz,be3.xyz),0.));
    /*if ( nrmbe1 > 1.E-6 ) trn4 = normalize(vec4(trn.xyz,0.))*length(be1);
    else trn4 = normalize(vec4(cross(be1.xyz, be3.xyz),0.));*/
    offset = vec4(trn4 * bary.w * thickness / 2.0f );
    gl_Position = bary-offset;
    color       = bcol;
    Nv          = be1.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();

    gl_Position =  bary+offset;
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

    gl_Position = bary- head_arrow_width * offset + 0.55 * be3;
    color       = bcol;
    Nv          = be1.xyz; //bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();

    gl_Position =  bary+ head_arrow_width * offset + 0.55 * be3;
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
    vec4 be2 = normalize(vec4(cross(be1.xyz,be3.xyz),0.));//length(be3)*vec4(trn.xyz,0.);
    vec3 nbe1 = normalize(be1.xyz).xyz;
    vec3 nbe2 = normalize(be2.xyz).xyz;
    be1 = normalize(be1);
    const float thickness = 0.015f;
    vec4 offset = vec4(be1 * bary.w * thickness / 2.0f );

    gl_Position = bary-offset;
    color       = bcol;
    Nv          = nbe2;//bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();

    gl_Position =  bary+offset;
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
    offset = vec4(be2 * bary.w * thickness / 2.0f );

    gl_Position = bary-offset;
    color       = bcol;
    Nv          = nbe1;//bnorm.xyz;
    P           = vbary.xyz;
    EmitVertex();

    gl_Position = bary+offset;
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

/* Calcul de la coordonee barycentrique u pour uniformisation des samples sur le maillage

    Reference : Efficient barycentric point sampling on meshes, J. Portsmouth, Pacific Graphics 2017
*/
float compute_u(float zeta_u, float zeta_v)
{
    const float tol=5.0E-3;
    float r = uniform_rand(); float l = (2.0f*zeta_u - zeta_v)/3.0;
    const int maxIter = 20; float u = 0.5; int n = 0;
    while ( n++ < maxIter )
    {
        float u1 = 1.0f-u;
        float P  = u*(2.0-u) - l*u*u1*u1 - r;
        float Pd = max(u1*(2.0f + l*(3.0f*u-1.0f)), 1.E-8);
        float du = max(min(P/Pd,0.25), -0.25); u -= du;
        u = max(min(u, 1.f-1.E-8),1.E-8);
        if ( abs(du) < tol ) break;
    }
    return u;
}

/* Calcul de la coordonee barycentrique v pour uniformisation des samples sur le maillage

    Reference : Efficient barycentric point sampling on meshes, J. Portsmouth, Pacific Graphics 2017
*/
float compute_v( float phi_u, float phi_v, float u )
{
    float r = uniform_rand();
    const float epsilon = 1.E-6;
    if ( abs(phi_v) < epsilon ) return (1.0f-u)*r;
    float tau = 1.f/3.f - (1.f + (u-1.f/3.f)*phi_u)/phi_v;
    float tmp = tau + u - 1.f;
    float q   = sqrt(tau*tau*(1.f-r) + tmp*tmp*r);
    return tau <= 0.5f*(1.f-u) ? tau + q : tau - q;
}
/*
    Genere un vecteur a l'interieur d'un triangle par inversion sampling ( par sommet barycentrique )

    Reference : Efficient barycentric point sampling on meshes, J. Portsmouth, Pacific Graphics 2017
*/
void generate_vector_inside_triangle()
{
    //float zeta_u = uniform_rand(), zeta_v = uniform_rand();
    float u = compute_u(1.f, 1.f);
    float v = compute_v(1.f, 1.f, u);
    float w = 1.f-u-v;
    vec4 world_coord = w*vertex[0].position + u * vertex[1].position + v*vertex[2].position;

    vec4 c = w*vertex[0].color+u*vertex[1].color+v*vertex[2].color;
    float f = length(c.rgb - vec3(0.5,0.5,0.5))*1.154700543;
    f = clamp(f, 0.0f, 1.0f);
    vec3 val = vec3(texture(colormap, f));
    c = vec4(val.r, val.g, val.b, 1.);

    vec4 nr = w*vertex[0].normal+u*vertex[1].normal+v*vertex[2].normal;
    vec4 be3= w*vertex[0].e3    +u*vertex[1].e3    +v*vertex[2].e3;

    if ( style_arrow == 0 ) draw_3d_arrow(world_coord, be3, nr, c);
    else if ( style_arrow == 1 ) draw_flat_arrow(world_coord, be3, nr, c);
    else if ( style_arrow == 2 ) draw_tetra_arrow(world_coord, be3, nr, c);
}

void generate_uniform_field()
{
	const float bary_ponderation = 1./3.;
	// Les trois sommets du triangle dans le repere World.
	vec4 TP1 = vertex[0].position;
    vec4 TP2 = vertex[1].position;
    vec4 TP3 = vertex[2].position;
    // Calcul du barycentre :
    vec4 Bary= bary_ponderation * (TP1 + TP2 + TP3);
    Bary.w = 1.;
    // Calcul de la surface du triangle :
    float surf_trig = 0.5*length(cross(TP2.xyz-TP1.xyz, TP3.xyz-TP1.xyz));
    // Calcul probabilite apparition d'un vecteur vitesse pour le triangle :
    float proba_trig = surf_trig*density;
    // Prendre une graine aleatoire en fonction du barycentre dans le repere monde pour garder 
    // une coherence visuelle lors des mouvements de la camera :
    seed = abs(int((length(Bary.xy)+length(Bary.yz)+length(Bary.xz))*16546.1415));
    // On se retrouve dans deux cas :
    //    1. Soit la probabilite est superieure a un : dans ce cas, on a au moins un vecteur dans le triangle, voire plus
    //       et on fait une boucle pour tirage aleatoire a l'interieur du triangle avec un nombre d'iteration egale
    //       au rapport de surface entre triangle et cellule cartesienne + 0.5...
    //    2. Soit la probabilite est inferieure a un. On fait un premier tirage pour savoir si un vecteur devra etre genere
    //       ou non a l'interieur de ce triangle :
    if ( proba_trig > 1. )
    {
        float uniform_value = uniform_rand();
        int niter = min(int(proba_trig+0.5),6);
        for ( int i = 0; i < niter; ++i )
            generate_vector_inside_triangle();
    }
    else
    {
        // Tirage probabilite qu'un point du triangle genere un vecteur :
        float uniform_value = rand(Bary.xyz); //uniform_rand();
        if ( uniform_value < proba_trig )
            generate_vector_inside_triangle();
    }
}

void main()
{
    if ( show_surface == 1 ) show_trig_surf();        // Emission du triangle de base

    // Emission des triangles pour le champs de vecteur
//    if ( density < 1.E-6 )
	if (density < 1.E-8)
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
