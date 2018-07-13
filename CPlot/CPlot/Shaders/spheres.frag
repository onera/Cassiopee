// Sphere billboarding
// Remplace un quad avec texcoord (0,0)-(1,1) par une sphere

varying vec4 initColor;
varying vec3 P;

void main()
{
    // r^2 = (x - x0)^2 + (y - y0)^2 + (z - z0)^2
    float x = gl_TexCoord[0].x-0.5;
    float y = gl_TexCoord[0].y-0.5;
    float zz = 0.25 - x*x - y*y;

    if (zz <= 0.) discard;
	
    float z = sqrt(zz);
    //float r = x*x+y*y;
    vec3 N = normalize(vec3(x, y, z));
    
    // light
    vec3 E = normalize(-P);
    vec3 L = normalize(gl_LightSource[0].position.xyz-P);
    float dotNL = dot(N, L);
    if (dotNL < 0.) { N = -N; dotNL = -dotNL; }

    vec3 R = normalize(-reflect(L,N));
    vec4 Iamb = gl_LightSource[0].ambient; 
    vec4 Idiff = gl_LightSource[0].diffuse*max(dotNL, 0.0);
    vec4 Ispec = gl_LightSource[0].specular * pow(max(dot(R,E),0.0),0.25 * gl_FrontMaterial.shininess);
    vec4 col = Iamb + initColor * Idiff + Ispec;

    // depth modification
    float depth = gl_FragCoord.z;
    depth += -z*0.0001;
	
    gl_FragDepth = depth;
    gl_FragColor = col;
    gl_FragColor.a = initColor.a;
    
}
