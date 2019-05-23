#version 150 compatibility

in vec4 color;
in vec3 Nv;
in vec3 P;
//in vec4 vert;

uniform int lightOn;
uniform float blend;
uniform int shadow;
uniform sampler2D ShadowMap;

void main()
{
  vec4 color2   = vec4(color.r, color.g, color.b, blend);
  vec4 color3   = vec4(0.25,0.25,0.25,blend);
  vec3 N = normalize(Nv);
  vec3 L = normalize(gl_LightSource[0].position.xyz-P);

  // light on!
  if (lightOn == 1)
  {
    vec3 E, R;
    vec4 Iamb, Idiff, Ispec;    
    E = normalize(-P);
    if (dot(N, L) < 0.) N = -N;

    R = normalize(-reflect(L,N));
    Iamb = gl_LightSource[0].ambient; 
    Idiff = gl_LightSource[0].diffuse*max(dot(N,L), 0.0);
    Ispec = gl_LightSource[0].specular * 0.1 * pow(max(dot(R,E),0.0),0.2 * gl_FrontMaterial.shininess);

    vec4 col = Iamb + color2 * Idiff + Ispec;
    vec4 colg = Iamb + color3 * Idiff + Ispec;
    color2 = clamp(col, 0., 1.);
    color3 = clamp(colg,0., 1.);
    color2.a = blend;
    color3.a = blend;
  }
  
  float shadowValue = 1.;
  if (shadow == 128.5)
  {
     // Coords -> texCoords
     vec4 ShadowCoord = gl_TextureMatrix[0] * color;
     vec4 shadowCoordinateW = ShadowCoord / ShadowCoord.w;

     // Used to lower moire pattern and self-shadowing
     float dotNL = dot(N, L);
     shadowCoordinateW.z -= (abs(dotNL)+0.1)*0.00001;

     // Z buffer du point dans la texture rendu du pt de vue de la lumiere
     float distanceFromLight = texture2D(ShadowMap, shadowCoordinateW.st).r;
     float s = shadowCoordinateW.s;
     float t = shadowCoordinateW.t;      
     if (ShadowCoord.w > 0.0 && s > 0.001 && s < 0.999 && t > 0.001 && t < 0.999)
       shadowValue = distanceFromLight < shadowCoordinateW.z ? 0.5 : 1.0;
  }
  if (color.a < 0.6) {
     gl_FragColor = shadowValue * color3;
     gl_FragColor.a = blend;
  } else {
    float b = ( color.a < blend ? color.a : blend );
    gl_FragColor = color2;
    gl_FragColor.a = b;
  }
}
