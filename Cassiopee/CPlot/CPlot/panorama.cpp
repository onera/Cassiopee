/*    
    Copyright 2013-2024 Onera.

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

# include "kcore.h"
# include "cplot.h"
# include "Data.h"

using namespace K_FLD;

//=============================================================================
/* 
   panorama: make a panorama from 6 cube images.
 */
//=============================================================================
PyObject* K_CPLOT::panorama(PyObject* self, PyObject* args)
{
  
  printf("Entering panorama\n");
  Data* d = Data::getInstance();

  // plutot le faire dans un dispay avec post processing

  // Create a texture to store the final 2:1 image
  E_Int height = 800;
  E_Int width = 2*height;
  d->createFrameBufferTexture();

  // Create and load the 6 cube map textures
  E_Int w, h;
  GLuint tex1, tex2, tex3, tex4, tex5, tex6;
  d->createImageTexture("front.png", tex1, w, h, false);
  d->createImageTexture("front.png", tex1, w, h, false);
  d->createImageTexture("front.png", tex1, w, h, false);
  d->createImageTexture("front.png", tex1, w, h, false);
  d->createImageTexture("front.png", tex1, w, h, false);

  // Load shader

  // Render and apply the shader

  // save the texture in a file

  return Py_BuildValue("i", KSUCCESS);
}
