/* XRay shader */ 
varying vec3 N;
varying vec3 I;
varying vec4 color;

uniform float EdgeFalloff;
uniform float intensity;

void main()
{
    float opacity = dot(normalize(N), normalize(I));
    opacity = abs(opacity);
    opacity = 1.0 - pow(opacity, EdgeFalloff);
    opacity = clamp(opacity, 0., 1.);
    
    //vec4 result = pow(color, vec4(1.0 / 2.2));
    gl_FragColor = color*(1.+opacity*2*(intensity-1));
    gl_FragColor.a = opacity;
}
