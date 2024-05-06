//
// Fragment shader for ray marching
//
varying vec4 color;
varying vec3 cam;
uniform sampler3D DensityTex;
uniform float Absorption;

void main()
{
    vec3 lightPos = vec3(cam.x, cam.y, cam.z-1.);
    vec3 lightIntensity = vec3(color)*12.;

    // diagonal of the cube
    const float maxDist = 1.73205080757;

    // Taille du voxel buffer
    const int numSamples = 96;
    const float scale = maxDist/float(numSamples);

    const int numLightSamples = 32;
    float lscale = maxDist / float(numLightSamples);
    float alscale = Absorption * lscale;
    float ascale = Absorption * scale;

    // assume all coordinates are in texture space
    vec3 pos = gl_TexCoord[0].xyz;
    vec3 eyeDir = normalize(pos-cam)*scale;

    // transmittance
    float T = 1.0;
    // in-scattered radiance
    vec3 Lo = vec3(0.0);

    for (int i=0; i < numSamples; i++)
    {
	// sample density
	float density = texture3D(DensityTex, pos).x;

	// skip empty space
	if (density > 0.0)
	{
	        // attenuate ray-throughput
		T *= 1.0-density*ascale;
		if (T <= 0.01) break;

		// point light dir in texture space
		vec3 lightDir = normalize(lightPos-pos)*lscale;

		// sample light
		float Tl = 1.0;	// transmittance along light ray
		vec3 lpos = pos + lightDir;

		for (int s=0; s < numLightSamples; ++s)
		{
		        float ld = texture3D(DensityTex, lpos).x;
			Tl *= 1.0-alscale*ld;

			if (Tl <= 0.01) break;

			lpos += lightDir;
		}

		vec3 Li = lightIntensity*Tl;

		Lo += Li*T*density*scale;
	}

	pos += eyeDir;
    }

    //gl_FragColor.r = gl_TexCoord[0].x;
    //gl_FragColor.g = gl_TexCoord[0].y;
    //gl_FragColor.b = 1; 

    gl_FragColor.rgb = Lo;
    gl_FragColor.a = 1.0-T;
    //gl_FragColor.a = color.a;
}
