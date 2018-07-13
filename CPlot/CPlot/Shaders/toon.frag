varying vec3 normal, lightDir;
varying vec4 color;

void main()
{
  float intensity;
  vec3 n;
  float factor;

  n = normalize(normal);
  intensity = dot(lightDir, n);

  if (intensity > 0.98) factor = 1.0;
  else if (intensity > 0.5) factor = 0.8;
  else if (intensity > 0.3) factor = 0.4;
  else factor = 0.0;
  gl_FragColor = color*vec4(factor, factor, factor, color.a);
}

