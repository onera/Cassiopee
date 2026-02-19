//
// Glass shader
//
const vec3 Xunitvec = vec3(1.0, 0.0, 0.0);
const vec3 Yunitvec = vec3(0.0, 1.0, 0.0);

varying vec4 color;
uniform float MixRatio2;
uniform float MixRatio;

// need to scale our framebuffer - it has a fixed width/height of 2048
uniform float FrameWidth;
uniform float FrameHeight;

uniform sampler2D EnvMap;
uniform sampler2D RefractionMap;

varying vec3 Normal;
varying vec3 EyeDir;
varying vec4 EyePos;

void main (void)
{
  float Depth = 0.2;
  vec3 BaseColor = vec3(color.r, color.g, color.b);

  vec3 LightPos = gl_LightSource[0].position.xyz;
  vec3 L = normalize(LightPos - EyeDir);
  float LightIntensity = dot(Normal, L);
  vec3 N = Normal;
  if (LightIntensity < 0.) 
  { N = -N; LightIntensity = -LightIntensity; }

  // Compute reflection vector
  vec3 reflectDir = reflect(EyeDir, N);

  // Compute altitude and azimuth angles
  vec2 index;
  index.y = dot(normalize(reflectDir), Yunitvec);
  reflectDir.y = 0.0;
  index.x = dot(normalize(reflectDir), Xunitvec) * 0.5;

  // Translate index values into proper range
  if (reflectDir.z >= 0.0)
    index = (index + 1.0) * 0.5;
  else
  {
    index.t = (index.t + 1.0) * 0.5;
    index.s = (-index.s) * 0.5 + 1.0;
  }
  // if reflectDir.z >= 0.0, s will go from 0.25 to 0.75
  // if reflectDir.z <  0.0, s will go from 0.75 to 1.25, and
  // that's OK, because we've set the texture to wrap.
  
  // Do a lookup into the environment map.
  vec3 envColor = vec3 (texture2D(EnvMap, index));
    
  // calc fresnels term.  This allows a view dependant blend of reflection/refraction
  float fresnel = abs(dot(normalize(EyeDir), N));
  fresnel *= MixRatio2;
  fresnel = clamp(fresnel, 0.1, 0.9);

  // calc refraction
  vec3 refractionDir = normalize(EyeDir) - normalize(N);

  // Scale the refraction so the z element is equal to depth
  float depthVal = Depth / -refractionDir.z;
	
  // perform the div by w
  float recipW = 1.0 / EyePos.w;
  vec2 eye = EyePos.xy * vec2(recipW);

  // calc the refraction lookup
  index.s = (eye.x + refractionDir.x * depthVal);
  index.t = (eye.y + refractionDir.y * depthVal);
	
  // scale and shift so we're in the range 0-1
  index.s = index.s / 2.0 + 0.5;
  index.t = index.t / 2.0 + 0.5;
	
  // as we're looking at the framebuffer, we want it clamping at the edge of the rendered scene, not the edge of the texture,
  // so we clamp before scaling to fit
  float recip1k = 1.0 / 2048.0;
  index.s = clamp(index.s, 0.0, 1.0 - recip1k);
  index.t = clamp(index.t, 0.0, 1.0 - recip1k);
	
  // scale the texture so we just see the rendered framebuffer
  index.s = index.s * FrameWidth * recip1k;
  index.t = index.t * FrameHeight * recip1k;
	
  vec3 RefractionColor = vec3(texture2D(RefractionMap, index));
    
  // Add lighting to base color and mix
  vec3 base = LightIntensity * BaseColor;
  envColor = mix(envColor, RefractionColor, fresnel);
  envColor = mix(envColor, base, MixRatio);
    
  //envColor = RefractionColor + 0.000001*envColor;

  gl_FragColor = vec4(envColor, color.a);
}