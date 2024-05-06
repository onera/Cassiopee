//
// Fragment shader for procedural bricks (with bumps)
//
uniform vec3 MortarColor;
uniform vec2 BrickSize;
uniform vec2 BrickPct;
uniform int shadow;
uniform sampler2D ShadowMap;

varying vec2 MCposition;
varying vec3 Nv;
varying vec3 P;
varying vec4 initColor;
varying vec4 vertex;

void main(void)
{
    vec3 color;
    vec2 position, useBrick;
    vec3 BrickColor = vec3(initColor.r, initColor.g, initColor.b);   
    position = MCposition / BrickSize;
    if (fract(position.y * 0.5) > 0.5) position.x += 0.5;
    position = fract(position);
    useBrick = step(position, BrickPct);
    color = mix(MortarColor, BrickColor, useBrick.x * useBrick.y);

     // Light
     vec3 N = normalize(Nv);
     vec3 E = normalize(-P);
     vec3 L = normalize(gl_LightSource[0].position.xyz-P);

     // bump
     float nx = 0., ny = 0.;
     if (position.x > 0.7) nx = (position.x-0.7);
     if (position.x < 0.2) nx = (position.x-0.2);
     if (position.y > 0.7) ny = (position.y-0.7);
     if (position.y < 0.2) ny = (position.y-0.2);
     N = N + 6.*gl_NormalMatrix * vec3(nx, ny ,0.);
     
     N = normalize(N);    
     if (dot(N, L) < 0.) N = -N;

     vec3 R = normalize(-reflect(L,N));
     vec4 Iamb = gl_LightSource[0].ambient;
     vec4 Idiff = gl_LightSource[0].diffuse*max(dot(N,L), 0.0);
     vec4 Ispec = gl_LightSource[0].specular * pow(max(dot(R,E),0.0), 0.2 * gl_FrontMaterial.shininess);
        
     vec4 col = vec4(color, initColor.a);   
     vec4 col2 = Iamb + col * Idiff + Ispec;
     col2 = clamp(col2, 0., 1.);

     float shadowValue = 1.;
     if (shadow > 0)
     {
     // Coords -> texCoords
     vec4 ShadowCoord = gl_TextureMatrix[0] * vertex;
     vec4 shadowCoordinateW = ShadowCoord / ShadowCoord.w;

     // Used to lower moire pattern and self-shadowing
     //shadowCoordinateW.z -= 0.00001;
     float dotNL = dot(N, L);
     shadowCoordinateW.z -= (abs(dotNL)+0.1)*0.00001;

     // Z buffer du point dans la texture rendu du pt de vue de la lumiere
     float distanceFromLight = texture2D(ShadowMap, shadowCoordinateW.st).r;
     float s = shadowCoordinateW.s;
     float t = shadowCoordinateW.t;      
     if (ShadowCoord.w > 0.0 && s > 0.001 && s < 0.999 && t > 0.001 && t < 0.999)
       shadowValue = distanceFromLight < shadowCoordinateW.z ? 0.5 : 1.0;
     }

     gl_FragColor = shadowValue * col2;
     gl_FragColor.a = initColor.a;
}