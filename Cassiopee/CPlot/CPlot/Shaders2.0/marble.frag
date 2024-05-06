//
// Marble shader
//
varying vec3  MCposition;
varying vec3 Nv;
varying vec3 P;
varying vec4 initColor;
varying vec4 vertex;

uniform sampler3D Noise;
uniform vec3 MarbleColor;
uniform vec3 VeinColor;
uniform int shadow;
uniform sampler2D ShadowMap;

// light is on?
uniform int lightOn;

void main()
{
    vec4 noisevec = texture3D(Noise, MCposition);

    float intensity = abs(noisevec[0] - 0.25) +
                      abs(noisevec[1] - 0.125) +
                      abs(noisevec[2] - 0.0625) +
                      abs(noisevec[3] - 0.03125);

    float sineval = sin(MCposition.y * 6.0 + intensity * 12.0) * 0.5 + 0.5;
    vec3 color = mix(VeinColor, MarbleColor, sineval);
    vec4 col;
    vec3 N = normalize(Nv);
    vec3 L = normalize(gl_LightSource[0].position.xyz-P);

    // Phong
    if (lightOn == 1)
    {
     vec3 E = normalize(-P);
     if (dot(N, L) < 0.) N = -N;
     vec3 R = normalize(-reflect(L,N));
     vec4 Iamb = gl_LightSource[0].ambient;
     vec4 Idiff = gl_LightSource[0].diffuse*max(dot(N,L), 0.0);
     vec4 Ispec = gl_LightSource[0].specular * pow(max(dot(R,E),0.0),0.2 * gl_FrontMaterial.shininess);
        
     col = vec4(color, initColor.a);   
     vec4 col2 =  Iamb + col * Idiff + Ispec;
     col = clamp(col2, 0., 1.);
     
    }
    else
    {
     col = vec4(color, initColor.a); 
    }

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

    gl_FragColor = shadowValue * col;
    gl_FragColor.a = initColor.a;
}