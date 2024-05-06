//
// Wood shader
//
varying vec3 MCposition;
varying vec3 Nv;
varying vec3 P;
varying vec4 initColor;
varying vec4 vertex;

uniform float bump;
uniform sampler3D Noise;
uniform int shadow;
uniform sampler2D ShadowMap;

uniform vec3  LightWoodColor;
uniform vec3  DarkWoodColor;
uniform float RingFreq;
uniform float LightGrains;
uniform float DarkGrains;
uniform float GrainThreshold;
uniform vec3  NoiseScale;
uniform float Noisiness;
uniform float GrainScale;

// light is on?
uniform int lightOn;

void main()
{
    vec3 noisevec = vec3(texture3D(Noise, MCposition*NoiseScale)*Noisiness);
    vec3 location = MCposition + noisevec;

    float dist = sqrt(location.x*location.x + location.z*location.z);
    dist *= RingFreq;
    float r = fract(dist + noisevec[0] + noisevec[1] + noisevec[2]) * 2.0;
    if (r > 1.0) r = 2.0 - r;
 
    vec3 color = mix(LightWoodColor, DarkWoodColor, r);

    r = fract((MCposition.x + MCposition.z) * GrainScale + 0.5);
    noisevec[2] *= r;
    if (r < GrainThreshold)
        color += LightWoodColor * LightGrains * noisevec[2];
    else
        color -= LightWoodColor * DarkGrains * noisevec[2];

    vec4 col;
    vec3 N = normalize(Nv);
    vec3 L = normalize(gl_LightSource[0].position.xyz-P);

    // Phong
    if (lightOn == 1)
    {
     vec3 E = normalize(-P);

     // bump
     noisevec = vec3(texture3D(Noise, MCposition));
     location = MCposition + noisevec;
     dist = sqrt(location.x*location.x + location.z*location.z);
     float r0 = fract(dist + 0.*noisevec[0] + 0.*noisevec[1] + 0.*noisevec[2]) * 2.0;
     if (r0 > 1.0) r0 = 2.0 - r0;
     //float r0 = (dist + 0.*noisevec[0] + 0.*noisevec[1] + 0.*noisevec[2])*0.8;
     //r0 = r0*r0-floor(r0*r0);

     noisevec = vec3(texture3D(Noise, MCposition+vec3(0.001,0.,0.)));
     location = MCposition+vec3(0.001,0.,0.) + noisevec;
     dist = sqrt(location.x*location.x + location.z*location.z);
     r = fract(dist + 0.*noisevec[0] + 0.*noisevec[1] + 0.*noisevec[2]) * 2.0;
     if (r > 1.0) r = 2.0 - r;
     //r = (dist + 0.*noisevec[0] + 0.*noisevec[1] + 0.*noisevec[2])*0.8;
     //r = r*r-floor(r*r);
     float nx = r - r0;

     noisevec = vec3(texture3D(Noise, MCposition+vec3(0.,0.001,0.)));
     location = MCposition+vec3(0.,0.001,0.) + noisevec;
     dist = sqrt(location.x*location.x + location.z*location.z);
     r = fract(dist + 0.*noisevec[0] + 0.*noisevec[1] + 0.*noisevec[2]) * 2.0;
     if (r > 1.0) r = 2.0 - r;
     //r = (dist + 0.*noisevec[0] + 0.*noisevec[1] + 0.*noisevec[2])*0.8;
     //r = r*r-floor(r*r);
     float ny = r - r0;

     noisevec = vec3(texture3D(Noise, MCposition+vec3(0.,0.,0.001)));
     location = MCposition+vec3(0.,0.,0.001) + noisevec;
     dist = sqrt(location.x*location.x + location.z*location.z);
     r = fract(dist + 0.*noisevec[0] + 0.*noisevec[1] + 0.*noisevec[2]) * 2.0;
     if (r > 1.0) r = 2.0 - r;
     //r = (dist + 0.*noisevec[0] + 0.*noisevec[1] + 0.*noisevec[2])*0.8;
     //r = r*r-floor(r*r);
     float nz = r - r0;

     N += bump*5.*gl_NormalMatrix * vec3(nx, ny ,nz);
     N = normalize(N);
     // end bump

     if (dot(N, L) < 0.) N = -N;
     vec3 R = normalize(-reflect(L,N));
     vec4 Iamb = gl_LightSource[0].ambient;
     vec4 Idiff = gl_LightSource[0].diffuse*max(dot(N,L), 0.0);
     vec4 Ispec = gl_LightSource[0].specular*pow(max(dot(R,E),0.0),0.1*gl_FrontMaterial.shininess);
        
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
