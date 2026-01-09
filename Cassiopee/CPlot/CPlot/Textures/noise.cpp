/*    
    Copyright 2013-2026 ONERA.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

// Coherent noise function over 1, 2 or 3 dimensions by K. Perlin
#include "../Data.h"
#include <cmath>
#include "kcore.h"

void init3DNoiseTexture(int texSize, GLubyte* texPtr);
void make3DNoiseTexture();

int Noise3DTexSize = 64;
GLubyte* Noise3DTexPtr;
//=============================================================================
void CreateNoise3D()
{
  make3DNoiseTexture();
  init3DNoiseTexture(Noise3DTexSize, Noise3DTexPtr);
}

//=============================================================================
void make3DNoiseTexture()
{
  int f, i, j, k, inc;
  int startFrequency = 4;
  int numOctaves = 4;
  double ni[3];
  double inci, incj, inck;
  int frequency = startFrequency;
  GLubyte* ptr;
  double amp = 0.5;
  K_NOISE::PDS data;
  
  Noise3DTexPtr = (GLubyte*)malloc(Noise3DTexSize*Noise3DTexSize*Noise3DTexSize*4*sizeof(GLubyte));

  for (f = 0, inc = 0; f < numOctaves; ++f, frequency *= 2, ++inc, amp *= 0.5)
  {
    K_NOISE::initPerlinNoise(frequency, data);
    ptr = Noise3DTexPtr;
    ni[0] = ni[1] = ni[2] = 0;
    
    inci = 1.0 / (Noise3DTexSize / frequency);
    for (i = 0; i < Noise3DTexSize; ++i, ni[0] += inci)
    {
      incj = 1.0 / (Noise3DTexSize / frequency);
      for (j = 0; j < Noise3DTexSize; ++j, ni[1] += incj)
      {
        inck = 1.0 / (Noise3DTexSize / frequency);
        for (k = 0; k < Noise3DTexSize; ++k, ni[2] += inck, ptr += 4)
          *(ptr + inc) = (GLubyte) (((K_NOISE::noise3(ni, data)+1.0)*amp)*128.0);
      }
    }
  }
}

//=============================================================================
void init3DNoiseTexture(int texSize, GLubyte* texPtr)
{
#ifdef __SHADERS__
  if (glewIsSupported("GL_EXT_texture3D") != 0)
  {
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, texSize, texSize, texSize, 
                 0, GL_RGBA, GL_UNSIGNED_BYTE, texPtr);
  }
#endif
  free(texPtr);
}
