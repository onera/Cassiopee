/*    
    Copyright 2013-2025 Onera.

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

#include "Images/png/png.h"
#include "../Data.h"

//=============================================================================
/*
  Create une texture a partir d'un fichier png.
  IN: filename: fichier png
  IN: mipmap: si true, cree une texture avec mipmaps (not used, always mipmap)
  OUT: tex: texture.
  OUT: width, height: nbre de pixels de la texture
*/
//=============================================================================
E_Int Data::createPngTexture(const char* filename, GLuint &tex, 
                             E_Int &width, E_Int &height, bool mipmap)
{
  Data* d = Data::getInstance();
  //GLint mipMap = GL_FALSE;
  //if (mipmap == true) mipMap = GL_TRUE;

  // Shader path
  char path[256*8];
  strcpy(path, d->ptrState->shaderPath);
#ifdef __SHADERS__
  strcat(path, filename);
#else
  strcpy(path, filename);
#endif

  // local path name
  char path2[256*8];
  //char* file = ptrState->file;
  char* lpn = ptrState->filePath;
  strcpy(path2, lpn);
  strcat(path2, "/");
  strcat(path2, filename);
  
  FILE* ptrFile = fopen(path, "rb"); // shader path
  if (!ptrFile) 
  { ptrFile = fopen(filename, "rb"); } // local cassiopee path
  if (!ptrFile)
  { ptrFile = fopen(path2, "rb"); } // loaded file path
  if (!ptrFile)
  { printf("Warning: CPlot: can't open texture file %s.\n", path2); 
    return 0; }
  
  png_structp png_ptr;
  png_byte color_type;
  png_byte bit_depth;

  // Header
  png_byte header[8];    // 8 is the maximum size that can be checked
  fread(header, 1, 8, ptrFile);
  if (png_sig_cmp(header, 0, 8))
  { printf("Warning: file %s is an invalid png file.\n", filename); 
    fclose(ptrFile); return 1; }
  
  /* initialize stuff */
  png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (!png_ptr)
  { printf("Warning: can not create png structures.\n"); fclose(ptrFile); return 1; }
  
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr)
  { printf("Warning: can not create png structures.\n"); fclose(ptrFile); return 1; }

  png_infop end_info_ptr = png_create_info_struct(png_ptr);
  if (!end_info_ptr)
  { printf("Warning: can not create png structures.\n"); fclose(ptrFile); return 1; }

  if (setjmp(png_jmpbuf(png_ptr)))
  { printf("Warning: can not create png structures.\n"); fclose(ptrFile); return 1; }
              
  png_init_io(png_ptr, ptrFile);
  png_set_sig_bytes(png_ptr, 8);

  png_read_info(png_ptr, info_ptr);

  width = png_get_image_width(png_ptr, info_ptr);
  height = png_get_image_height(png_ptr, info_ptr);
  color_type = png_get_color_type(png_ptr, info_ptr);
  bit_depth = png_get_bit_depth(png_ptr, info_ptr);
  //printf("color type %d\n", color_type);
  //printf("bit depth %d\n", bit_depth);
  
  if (bit_depth == 16) png_set_strip_16(png_ptr);
  if (color_type == PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png_ptr);

  // PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth.
  if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) png_set_expand_gray_1_2_4_to_8(png_ptr);

  if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS)) png_set_tRNS_to_alpha(png_ptr);
  // These color_type don't have an alpha channel then fill it with 0xff.
  if (color_type == PNG_COLOR_TYPE_RGB ||
      color_type == PNG_COLOR_TYPE_GRAY ||
      color_type == PNG_COLOR_TYPE_PALETTE)
    png_set_filler(png_ptr, 0xFF, PNG_FILLER_AFTER);

  if (color_type == PNG_COLOR_TYPE_GRAY ||
      color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
    png_set_gray_to_rgb(png_ptr);

  //int number_of_passes = png_set_interlace_handling(png_ptr);
  png_read_update_info(png_ptr, info_ptr);

  /* read file */
  if (setjmp(png_jmpbuf(png_ptr)))
  { printf("Warning: error during read of png image.\n"); fclose(ptrFile); return 1; }
  
  bool alpha; int components;
  E_Int size = png_get_rowbytes(png_ptr, info_ptr);
  if (size == width*3) { alpha = false; components = 3; }
  else { alpha = true; components = 4; }

  png_bytepp rows = new png_bytep[height];
  png_bytep image = new png_byte[height * width * components];
  
  for (E_Int i = 0; i < height; i++)
    rows[i] = image + (height - i - 1) * width * components;
  
  png_set_rows(png_ptr, info_ptr, rows);
  png_read_image(png_ptr, rows);

  //printf("reading file %s\n", filename);

  if (tex == 0)
  {
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
    glTexImage2D(GL_TEXTURE_2D, 0, components, width, height, 0, 
                 alpha ? GL_RGBA : GL_RGB, GL_UNSIGNED_BYTE, image);
  }
  else
  {
    glDeleteTextures(1, &tex);
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    //glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height,
    //                alpha ? GL_RGBA : GL_RGB, GL_UNSIGNED_BYTE, image);
    glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
    glTexImage2D(GL_TEXTURE_2D, 0, components, width, height, 0, 
    alpha ? GL_RGBA : GL_RGB, GL_UNSIGNED_BYTE, image);
  }

  png_read_end(png_ptr, end_info_ptr);
  png_destroy_read_struct(&png_ptr, &info_ptr, &end_info_ptr);
  delete [] image;
  delete [] rows;
  fclose(ptrFile);
  return 1;
}
