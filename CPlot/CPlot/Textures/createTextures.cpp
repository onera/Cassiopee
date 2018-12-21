/*    
    Copyright 2013-2019 Onera.

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

#include "../Data.h"
#ifdef hypot
#undef hypot
#endif

void CreateNoise3D();

//=============================================================================
/*
  Create texture for nodes. Cette texture node est utilise dans displayNode.
  Create a mipmap of textures.
*/
//=============================================================================
int Data::createNodeTexture(void)
{
  int texWidth = 256;
  int texHeight = 256;
  GLubyte *texPixels, *p;
  int texSize;
  int i, j;
  
  glGenTextures(1, &_texNodes);
  glBindTexture(GL_TEXTURE_2D, _texNodes);
   
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, 
                  GL_LINEAR_MIPMAP_NEAREST);

  texSize = texWidth*texHeight*4*sizeof(GLubyte);
  texPixels = (GLubyte *)malloc(texSize);
  if (texPixels == NULL)
  {
    _texNodes = 0;
    return 0;
  }

  p = texPixels;
  for (i = 0; i < texHeight; i++) 
  {
    for (j = 0; j < texWidth; j++) 
    {
      float dist =
        hypot((float)(i - (texHeight / 2.)),
              (float)(j - (texWidth / 2.)));
      
      int color =  (int)( 255.- (dist*1.8) );
      if (color < 0) color = 0;
      p[0] = color; p[1] = color; p[2] = color; p[3] = color;
      p += 4;
    }
  }

  gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, texWidth, texHeight, GL_RGBA, 
                    GL_UNSIGNED_BYTE, texPixels);
  
  free(texPixels);
  return 1;
}

//=============================================================================
/*
  Create texture for Noise3D.
  Noise3D est une texture 3D qui sert comme tirage aleatoire dans les shaders.
*/
//=============================================================================
int Data::createNoise3DTexture(void)
{
#ifdef __SHADERS__
  if (glewIsSupported("GL_EXT_texture3D") != 0)
  {
    glGenTextures(1, &_texNoise);
    glBindTexture(GL_TEXTURE_3D, _texNoise);
    CreateNoise3D();
  }
#endif
  return 1;
}

//=============================================================================
/*
  Create Frame Buffer texture.
  Cette texture sert a stocker une copie de l'ecran (courant)
*/
//=============================================================================
int Data::createFrameBufferTexture(void)
{
  glGenTextures(1, &_texFrameBuffer);
  glBindTexture(GL_TEXTURE_2D, _texFrameBuffer);
  glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, 
                   _frameBufferSize, _frameBufferSize, 0);
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  return 1;
}
