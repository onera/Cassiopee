#version 130
// Depth of field shader (post shader)
uniform sampler2D FrameBuffer;
uniform sampler2D depthMap;

uniform float focalDepth; // DOF: position de la focale
uniform float radius; // DOF: taille du rayon de blur 
uniform float ext; // DOF: extension du blur

uniform int toneMapping; // TONE: type de tone mapping
uniform float gamma; // TONE: gamma correction

uniform float sobelThreshold; // SOBEL: threshold for sobel

uniform float sharpenCoeff; // SHARPEN: coeff

uniform float ssaoRadius; // SSAO: radius for evaluation

// internes
vec2 poisson0, poisson1, poisson2, poisson3, poisson4;
vec2 poisson5, poisson6, poisson7;
vec2 maxCoC = vec2(radius, 2.*radius);
vec2 pixelSizeHigh;
vec4 sobelColor = vec4(0,0,0,1);

vec4 dof(vec2 coords)
{
  // pas de blur au point focal
  // puis lineaire avant et apres sur dist en depth
  float dist = 0.02/(2.*radius);
  vec4 finalColor;
  float discRadius, discRadiusLow, centerDepth;

  centerDepth = texture2D(depthMap, coords).r;
  centerDepth = (centerDepth - focalDepth)/dist;
  centerDepth = pow(centerDepth, ext+0.0001);
  centerDepth = clamp(centerDepth,-1.,1.);
  centerDepth = centerDepth*0.5 + 0.5;

  discRadius = abs(centerDepth*maxCoC.y - maxCoC.x);
  discRadiusLow = discRadius * 0.5;

  finalColor = vec4(0.0,0.0,0.0,0.0);
  float finalBlur = 0.;
  float blur, tapBlur;
    
  // -- 
  vec2 coordHigh = coords + (pixelSizeHigh*poisson0*discRadius);
  vec4 tapHigh = texture2D(FrameBuffer, coordHigh);
  float blurHigh = texture2D(depthMap, coordHigh).r;
  blurHigh = (blurHigh - focalDepth)/dist;
  blurHigh = clamp(blurHigh,-1.,1.);
  blurHigh = blurHigh*0.5 + 0.5;
  tapBlur = abs(blurHigh*2.0 - 1.0);
  blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
  finalColor += tapHigh * blur;
  finalBlur += blur;
      
  // -- 
  coordHigh = coords + (pixelSizeHigh*poisson1*discRadius);
  tapHigh = texture2D(FrameBuffer, coordHigh);
  blurHigh = texture2D(depthMap, coordHigh).r;
  blurHigh = (blurHigh - focalDepth)/dist;
  blurHigh = clamp(blurHigh,-1.,1.);
  blurHigh = blurHigh*0.5 + 0.5;
  tapBlur = abs(blurHigh*2.0 - 1.0);
  blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
  finalColor += tapHigh * blur;
  finalBlur += blur;

  // -- 
  coordHigh = coords + (pixelSizeHigh*poisson2*discRadius);
  tapHigh = texture2D(FrameBuffer, coordHigh);
  blurHigh = texture2D(depthMap, coordHigh).r;
  blurHigh = (blurHigh - focalDepth)/dist;
  blurHigh = clamp(blurHigh,-1.,1.);
  blurHigh = blurHigh*0.5 + 0.5;
  tapBlur = abs(blurHigh*2.0 - 1.0);
  blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
  finalColor += tapHigh * blur;
  finalBlur += blur;

  // -- 
  coordHigh = coords + (pixelSizeHigh*poisson3*discRadius);
  tapHigh = texture2D(FrameBuffer, coordHigh);
  blurHigh = texture2D(depthMap, coordHigh).r;
  blurHigh = (blurHigh - focalDepth)/dist;
  blurHigh = clamp(blurHigh,-1.,1.);
  blurHigh = blurHigh*0.5 + 0.5;
  tapBlur = abs(blurHigh*2.0 - 1.0);
  blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
  finalColor += tapHigh * blur;
  finalBlur += blur;

  // -- 
  coordHigh = coords + (pixelSizeHigh*poisson4*discRadius);
  tapHigh = texture2D(FrameBuffer, coordHigh);
  blurHigh = texture2D(depthMap, coordHigh).r;
  blurHigh = (blurHigh - focalDepth)/dist;
  blurHigh = clamp(blurHigh,-1.,1.);
  blurHigh = blurHigh*0.5 + 0.5;
  tapBlur = abs(blurHigh*2.0 - 1.0);
  blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
  finalColor += tapHigh * blur;
  finalBlur += blur;

  // -- 
  coordHigh = coords + (pixelSizeHigh*poisson5*discRadius);
  tapHigh = texture2D(FrameBuffer, coordHigh);
  blurHigh = texture2D(depthMap, coordHigh).r;
  blurHigh = (blurHigh - focalDepth)/dist;
  blurHigh = clamp(blurHigh,-1.,1.);
  blurHigh = blurHigh*0.5 + 0.5;
  tapBlur = abs(blurHigh*2.0 - 1.0);
  blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
  finalColor += tapHigh * blur;
  finalBlur += blur;

  // -- 
  coordHigh = coords + (pixelSizeHigh*poisson6*discRadius);
  tapHigh = texture2D(FrameBuffer, coordHigh);
  blurHigh = texture2D(depthMap, coordHigh).r;
  blurHigh = (blurHigh - focalDepth)/dist;
  blurHigh = clamp(blurHigh,-1.,1.);
  blurHigh = blurHigh*0.5 + 0.5;
  tapBlur = abs(blurHigh*2.0 - 1.0);
  blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
  finalColor += tapHigh * blur;
  finalBlur += blur;

  // -- 
  coordHigh = coords + (pixelSizeHigh*poisson7*discRadius);
  tapHigh = texture2D(FrameBuffer, coordHigh);
  blurHigh = texture2D(depthMap, coordHigh).r;
  blurHigh = (blurHigh - focalDepth)/dist;
  blurHigh = clamp(blurHigh,-1.,1.);
  blurHigh = blurHigh*0.5 + 0.5;
  tapBlur = abs(blurHigh*2.0 - 1.0);
  blur = ((blurHigh >= centerDepth) ? 1.0 : tapBlur);
  finalColor += tapHigh * blur;
  finalBlur += blur;

  return finalColor/finalBlur;
}

vec4 sobel(vec2 coords)
{
  vec4 he = vec4(0.);
  float dx = pixelSizeHigh.x;
  float dy = pixelSizeHigh.y;
  he -= texture2D(FrameBuffer, vec2(coords.x - dx, coords.y - dy));
  he -= texture2D(FrameBuffer, vec2(coords.x - dx, coords.y     ))*2.;
  he -= texture2D(FrameBuffer, vec2(coords.x - dx, coords.y + dy));
  he += texture2D(FrameBuffer, vec2(coords.x + dx, coords.y - dy));
  he += texture2D(FrameBuffer, vec2(coords.x + dx, coords.y     ))*2.;
  he += texture2D(FrameBuffer, vec2(coords.x + dx, coords.y + dy));

  vec4 ve = vec4(0.);
  ve -= texture2D(FrameBuffer, vec2(coords.x - dx, coords.y - dy));
  ve -= texture2D(FrameBuffer, vec2(coords.x     , coords.y - dy))*2.;
  ve -= texture2D(FrameBuffer, vec2(coords.x + dx, coords.y - dy));
  ve += texture2D(FrameBuffer, vec2(coords.x - dx, coords.y + dy));
  ve += texture2D(FrameBuffer, vec2(coords.x     , coords.y + dy))*2.;
  ve += texture2D(FrameBuffer, vec2(coords.x + dx, coords.y + dy));

  vec3 e = sqrt(he.rgb*he.rgb + ve.rgb*ve.rgb);

  return vec4(e, 1.);
}

// Narkowicz 2015, "ACES Filmic Tone Mapping Curve"
vec3 aces(vec3 x) 
{
  const float a = 2.51;
  const float b = 0.03;
  const float c = 2.43;
  const float d = 0.59;
  const float e = 0.14;
  return clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0, 1.0);
}

// Filmic Tonemapping Operators http://filmicworlds.com/blog/filmic-tonemapping-operators/
vec3 filmic(vec3 x) 
{
  vec3 X = max(vec3(0.0), x - 0.004);
  vec3 result = (X * (6.2 * X + 0.5)) / (X * (6.2 * X + 1.7) + 0.06);
  return pow(result, vec3(2.2));
}

// Uchimura 2017, "HDR theory and practice"
vec3 uchimura(vec3 x, float P, float a, float m, float l, float c, float b) 
{
  float l0 = ((P - m) * l) / a;
  float L0 = m - m / a;
  float L1 = m + (1.0 - m) / a;
  float S0 = m + l0;
  float S1 = m + a * l0;
  float C2 = (a * P) / (P - S1);
  float CP = -C2 / P;

  vec3 w0 = vec3(1.0 - smoothstep(0.0, m, x));
  vec3 w2 = vec3(step(m + l0, x));
  vec3 w1 = vec3(1.0 - w0 - w2);

  vec3 T = vec3(m * pow(x / m, vec3(c)) + b);
  vec3 S = vec3(P - (P - S1) * exp(CP * (x - S0)));
  vec3 L = vec3(m + a * (x - m));

  return T * w0 + L * w1 + S * w2;
}

vec3 uchimura(vec3 x) 
{
  const float P = 1.0;  // max display brightness
  const float a = 1.0;  // contrast
  const float m = 0.22; // linear section start
  const float l = 0.4;  // linear section length
  const float c = 1.33; // black
  const float b = 0.0;  // pedestal

  return uchimura(x, P, a, m, l, c, b);
}

// on doit compter les voisins dont le depth est plus petit
float ssao(vec2 coords)
{
  float dist = 0.02/(2.*radius);
  float discRadius, centerDepth;
  float occlusion = 0.;
  float cdepth = texture2D(depthMap, coords).r; // current fragment depth
  centerDepth = cdepth;
  //centerDepth = (centerDepth - focalDepth)/dist;
  //centerDepth = pow(centerDepth, ext+0.0001);
  //centerDepth = clamp(centerDepth,-1.,1.);
  //centerDepth = centerDepth*0.5 + 0.5;
  discRadius = radius;
  
  vec2 coordHigh = coords + (pixelSizeHigh*poisson0*discRadius);
  float depth = texture2D(depthMap, coordHigh).r;
  if (depth < cdepth) occlusion += 1.;

  coordHigh = coords + (pixelSizeHigh*poisson1*discRadius);
  depth = texture2D(depthMap, coordHigh).r;
  if (depth < cdepth) occlusion += 1.;

  coordHigh = coords + (pixelSizeHigh*poisson2*discRadius);
  depth = texture2D(depthMap, coordHigh).r;
  if (depth < cdepth) occlusion += 1.;

  coordHigh = coords + (pixelSizeHigh*poisson3*discRadius);
  depth = texture2D(depthMap, coordHigh).r;
  if (depth < cdepth) occlusion += 1.;

  coordHigh = coords + (pixelSizeHigh*poisson4*discRadius);
  depth = texture2D(depthMap, coordHigh).r;
  if (depth < cdepth) occlusion += 1.;

  coordHigh = coords + (pixelSizeHigh*poisson5*discRadius);
  depth = texture2D(depthMap, coordHigh).r;
  if (depth < cdepth) occlusion += 1.;

  coordHigh = coords + (pixelSizeHigh*poisson6*discRadius);
  depth = texture2D(depthMap, coordHigh).r;
  if (depth < cdepth) occlusion += 1.;

  coordHigh = coords + (pixelSizeHigh*poisson7*discRadius);
  depth = texture2D(depthMap, coordHigh).r;
  if (depth < cdepth) occlusion += 1.;

  occlusion = occlusion / 8.;
  return occlusion;
}

// sharpen
vec4 sharpen(vec2 coords, float coeff) 
{
  float dx = pixelSizeHigh.x;
  float dy = pixelSizeHigh.y;
  vec4 sum = vec4(0.0);
  sum += -coeff * texture2D(FrameBuffer, coords + vec2( -1.0 * dx , 0.0 * dy));
  sum += -coeff * texture2D(FrameBuffer, coords + vec2( 0.0 * dx , -1.0 * dy));
  sum += (1.+4.*coeff) * texture2D(FrameBuffer, coords + vec2( 0.0 * dx , 0.0 * dy));
  sum += -coeff * texture2D(FrameBuffer, coords + vec2( 0.0 * dx , 1.0 * dy));
  sum += -coeff * texture2D(FrameBuffer, coords + vec2( 1.0 * dx , 0.0 * dy));
  return sum;
}

//=====
void main()
{
  pixelSizeHigh = 1.0/textureSize(FrameBuffer, 0);
 
  // poisson-distributed positions
  poisson0 = vec2( 0.0, 0.0);
  poisson1 = vec2( 0.527837,-0.08586);
  poisson2 = vec2(-0.040088, 0.536087);
  poisson3 = vec2(-0.670445,-0.179949);
  poisson4 = vec2(-0.419418,-0.616039);
  poisson5 = vec2( 0.440453,-0.639399);
  poisson6 = vec2(-0.757088, 0.349334);
  poisson7 = vec2( 0.574619, 0.685879);

  // color
  vec4 color1, color2, color;
  
  // dof pass
  if (radius > 0.) color1 = dof(gl_TexCoord[0].xy);
  else color1 = texture2D(FrameBuffer, gl_TexCoord[0].xy);

  // sobel pass
  if (sobelThreshold > 0.) color2 = sobel(gl_TexCoord[0].xy);
  else color2 = vec4(0.);

  // mix colors
  if (sobelThreshold > 0.)
  {
    float ct = max(color2.r, color2.g);
    ct = max(ct, color2.b);
    //ct = smoothstep(0.4, 0.8, ct);
    if (ct > sobelThreshold) color = sobelColor;
    else 
    color = mix(color1, 1.-color2, 0.3); // pastel colors
    //color = mix(color1,sobelColor,(1./sobelThreshold)*ct+1.);
  }
  else color = color1;

  vec4 res = color;
  
  // sharpen
  if (sharpenCoeff > 0.)
  {
    res = sharpen(gl_TexCoord[0].xy, sharpenCoeff);
    res.a = color.a;
  }
  
  // ssao
  if (ssaoRadius > 0.)
  {
    float occlusion = ssao(gl_TexCoord[0].xy);
    res.rgb = (1.-occlusion)*res.rgb;
  }

  // tone pass
  if (toneMapping == 0) // pure gamma correction
  {
    res = pow(res, vec4(1.0 / gamma));
    res.a = color.a;
  }
  else if (toneMapping == 1) // ACES
  {
    res.rgb = aces(res.rgb);
    res.rgb = pow(res.rgb, vec3(1.0 / gamma));
    res.a = color.a;
  }
  else if (toneMapping == 2) // Filmic
  {
    res.rgb = filmic(res.rgb);
    res.rgb = pow(res.rgb, vec3(1.0 / gamma));
    res.a = color.a;
  }
  else if (toneMapping == 3) // Uchimura
  {
    res.rgb = uchimura(res.rgb);
    res.rgb = pow(res.rgb, vec3(1.0 / gamma));
    res.a = color.a;
  }
  //res.rgb = vec3(1.-occlusion);
  
  gl_FragColor = res;
}
