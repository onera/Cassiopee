#version 150 compatibility
layout(triangles) in;
layout(line_strip, max_vertices=100) out;

in Vertex {
    vec4 P0;
    vec4 P1;
    vec4 normal;
    vec4 translation;
    vec4 color;
} vertex[];

out vec4 color;
out float gAlpha;

uniform sampler1D colormap;
uniform float density;

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

void simple_line_draw()
{
    int i;
    float ust = 0.33333333333;
    vec4 bary1 = gl_ModelViewProjectionMatrix * (ust*(vertex[0].P0+vertex[1].P0+vertex[2].P0));
    vec4 bary2 = gl_ModelViewProjectionMatrix * (ust*(vertex[0].P1+vertex[1].P1+vertex[2].P1));
    vec4 bcol = ust*(vertex[0].color+vertex[1].color+vertex[2].color);
    float f = length(bcol.rgb - vec3(0.5,0.5,0.5))*1.154700543;
    f = clamp(f, 0.0f, 1.0f);
    vec3 val = vec3(texture(colormap, f));
    bcol = vec4(val.r, val.g, val.b, 1.);
    
    gl_Position = bary1;
    color  = bcol;
    gAlpha = 0.8f;
    EmitVertex();

    gl_Position = bary2;
    color  = bcol;
    gAlpha = 0.f;
    EmitVertex();
    EndPrimitive();    
}

/* Calcul de la coordonee barycentrique u pour uniformisation des samples sur le maillage

    Reference : Efficient barycentric point sampling on meshes, J. Portsmouth, Pacific Graphics 2017
*/
float compute_u(float zeta_u, float zeta_v, const float tol = 5.0E-3 )
{
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
void generate_line_inside_triangle()
{
    //float zeta_u = uniform_rand(), zeta_v = uniform_rand();
    float u = compute_u(1.f, 1.f);
    float v = compute_v(1.f, 1.f, u);
    float w = 1.f-u-v;
    vec4  wp= w*vertex[0].P0 + u * vertex[1].P0 + v*vertex[2].P0;
    vec4  p= gl_ModelViewProjectionMatrix * wp;

    vec4 c = w*vertex[0].color+u*vertex[1].color+v*vertex[2].color;
    float f = length(c.rgb - vec3(0.5,0.5,0.5))*1.154700543;
    f = clamp(f, 0.0f, 1.0f);
    vec3 val = vec3(texture(colormap, f));
    c = vec4(val.r, val.g, val.b, 1.);

    vec4 t = w*vertex[0].translation+u*vertex[1].translation+v*vertex[2].translation;

    gl_Position = p;
    color = c;
    gAlpha = 0.8f;
    EmitVertex();

    gl_Position = p+t;
    color = c;
    gAlpha = 0.f;
    EmitVertex();
    EndPrimitive();
}

void uniform_line_draw()
{
    const float bary_ponderation = 1./3.;
    // Les trois sommets du triangle dans le repere World.
    vec4 TP1 = vertex[0].P0;
    vec4 TP2 = vertex[1].P0;
    vec4 TP3 = vertex[2].P0;
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
        int niter = min(int(proba_trig+0.5),50);
        for ( int i = 0; i < niter; ++i )
            generate_line_inside_triangle();
    }
    else
    {
        // Tirage probabilite qu'un point du triangle genere un vecteur :
        float uniform_value = rand(Bary.xyz); //uniform_rand();
        if ( uniform_value < proba_trig )
            generate_line_inside_triangle();
    }
/*    int i;

    vec4 v0P0 = gl_ModelViewProjectionMatrix * vertex[0].P0;
    vec4 v1P0 = gl_ModelViewProjectionMatrix * vertex[1].P0;
    vec4 v2P0 = gl_ModelViewProjectionMatrix * vertex[2].P0;
    vec4 e1 = v1P0 - v0P0;
    vec4 e2 = v2P0 - v0P0;
    vec3 unrm = cross(e1.xyz,e2.xyz);
    float nrm = length(unrm);
    //float nrm = sqrt(dot(unrm.xyz,unrm.xyz));
    //float nrm_e1 = sqrt(e1.x*e1.x+e1.y*e1.y+e1.z*e1.z);
    //float nrm_e2 = sqrt(e2.x*e2.x+e2.y*e2.y+e2.z*e2.z);
    //int ni = int(nrm_e1*density);
    //int nj = int(nrm_e2*density);
    float proba = nrm*density;
    int n = int(0.5f*(-1.f+sqrt(1.+8.f*proba)));
    n = ( n > 4 ? 4 : n );
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
            vec4 t = psi*vertex[0].translation+ki*vertex[1].translation+te*vertex[2].translation;

            gl_Position = p;
            color = c;
            gAlpha = 0.8f;
            EmitVertex();

            gl_Position = p+t;
            color = c;
            gAlpha = 0.f;
            EmitVertex();
            EndPrimitive();
        }
    } }
    else {
    //if ( n == 0 ) {
        float ust = (1./3.);
        // On va tirer une proba pour voir si on emet un vecteur :
        vec4 bary1 = ust*(v0P0+v1P0+v2P0);
        //float proba = (nrm_e1+nrm_e2)*0.5f*density;
        float r = 0.5*rand(bary1.xy);
        if (r < proba)
        {
            vec4 t = ust*(vertex[0].translation+vertex[1].translation+vertex[2].translation);
            //vec4 bvert= ust*(vertex[0].position+vertex[1].position+vertex[2].position);
            vec4 bcol = ust*(vertex[0].color+vertex[1].color+vertex[2].color);
            float f = length(bcol.rgb - vec3(0.5,0.5,0.5))*1.154700543;
            f = clamp(f, 0.0f, 1.0f);
            vec3 val = vec3(texture(colormap, f));
            bcol = vec4(val.r, val.g, val.b, 1.);

            gl_Position = bary1;
            color  = bcol;
            gAlpha = 0.8f;
            EmitVertex();

            gl_Position =  bary1+t;
            color  = bcol;
            gAlpha = 0.f;
            EmitVertex();
            EndPrimitive();
        }
    }
*/
}

void main()
{
    if ( density < 1.E-8 ) simple_line_draw();
    else uniform_line_draw();
}
