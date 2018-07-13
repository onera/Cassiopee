/*
    Phong (two sides) + shadow
*/
varying vec3 Nv;
varying vec3 P;
varying vec4 color;
varying vec4 vertex;
uniform float specularFactor;
uniform int shadow;
uniform sampler2D ShadowMap;

void main (void)
{ 
  vec3 N = normalize(Nv);
  vec3 E = normalize(-P);
  vec3 L = normalize(gl_LightSource[0].position.xyz-P);
  float dotNL = dot(N, L);
  if (dotNL < 0.) { N = -N; dotNL = -dotNL; }

  vec3 R = normalize(-reflect(L,N));
  vec4 Iamb = gl_LightSource[0].ambient; 
  vec4 Idiff = gl_LightSource[0].diffuse*max(dotNL, 0.0);
  vec4 Ispec = (specularFactor*specularFactor)*gl_LightSource[0].specular*pow(max(dot(R,E),0.0),0.25*gl_FrontMaterial.shininess);
  vec4 col = Iamb + color*Idiff + Ispec;
  col = clamp(col, 0., 1.);
  float shadowValue = 1.;
        
  if (shadow > 0)
  {
  // Coords -> texCoords
  vec4 ShadowCoord = gl_TextureMatrix[0] * vertex;
  vec4 shadowCoordinateW = ShadowCoord / ShadowCoord.w;

  // Used to lower moire pattern and self-shadowing
  //shadowCoordinateW.z -= 0.00001;
  shadowCoordinateW.z -= (abs(dotNL)+0.1)*0.00001;

  // Z buffer du point dans la texture rendu du pt de vue de la lumiere
  float distanceFromLight = texture2D(ShadowMap, shadowCoordinateW.st).r;
  float s = shadowCoordinateW.s;
  float t = shadowCoordinateW.t;      
  if (ShadowCoord.w > 0.0 && s > 0.001 && s < 0.999 && t > 0.001 && t < 0.999)
      shadowValue = distanceFromLight < shadowCoordinateW.z ? 0.5 : 1.0;
  }

  gl_FragColor = shadowValue * col;
  gl_FragColor.a = color.a;
}
