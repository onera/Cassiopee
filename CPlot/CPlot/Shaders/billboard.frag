// Texture billboarding
// avec eclairage spherique

uniform sampler2D BTex;
uniform int Ni;
uniform int Nj;

varying vec4 initColor;
varying vec3 P;

void main()
{
    vec4 color;
    vec3 c1; vec3 c2;
    vec2 index = gl_TexCoord[0].st;

    // r^2 = (x - x0)^2 + (y - y0)^2 + (z - z0)^2
    float x = index.x-0.5;
    float y = index.y-0.5;
    float zz = 0.25 - x*x - y*y;
        
    // Choix de la texture base sur r
    int e1, e2;
    int part = int(initColor.r*65280.);
    e2 = part/255;
    e1 = part-255*e2;
    float r = float(e1)/255.;
    float cr = float(e2)/255.;

    float di = 1./float(Ni);
    float dj = 1./float(Nj);
    int p = int(r*float(Ni)*float(Nj));
    int pj = p/Ni;
    int pi = p - pj*Ni;
    index.x = index.x * di + float(pi)*di;
    index.y = index.y * dj + float(pj)*dj;

    color = texture2D(BTex, index);
    c1 = vec3(color);
    c2 = vec3(initColor); c2.r = cr; // forced
    c1 = c1 * c2;
    color = vec4(c1, color.a);
    //color = vec4(r,r,r,color.a);
    //if (color.a < 0.000001) discard;

    // eclairage
    vec4 col;

    // eclairage spherique supprime:
    zz = -1.;

    if (zz <= 0.) col = color;
    else
    {
        float z = sqrt(zz);
        vec3 N = normalize(vec3(x, y, z));
        vec3 E = normalize(-P);
        vec3 L = normalize(gl_LightSource[0].position.xyz-P);
        float dotNL = dot(N, L);
        if (dotNL < 0.) { N = -N; dotNL = -dotNL; }

        vec3 R = normalize(-reflect(L,N));
        vec4 Iamb = gl_LightSource[0].ambient; 
        vec4 Idiff = gl_LightSource[0].diffuse*max(dotNL, 0.0);
        //vec4 Ispec = gl_LightSource[0].specular * pow(max(dot(R,E),0.0),0.25 * gl_FrontMaterial.shininess);
        col = Iamb + color * Idiff;
    }

    col.a = clamp(color.a - (1.-initColor.a), 0., 1.);
    gl_FragColor = col;
}
