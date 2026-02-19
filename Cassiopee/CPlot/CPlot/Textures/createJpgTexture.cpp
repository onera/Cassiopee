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

#include "../Data.h"

#include "Images/libjpeg/jpeglib.h"

struct my_error_mgr 
{
  struct jpeg_error_mgr pub;	/* "public" fields */
  //jmp_buf setjmp_buffer;	/* for return to caller */
};
typedef struct my_error_mgr* my_error_ptr;
void my_error_exit (j_common_ptr cinfo)
{
  /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
  //my_error_ptr myerr = (my_error_ptr) cinfo->err;

  /* Always display the message. */
  /* We could postpone this until after returning, if we chose. */
  (*cinfo->err->output_message) (cinfo);
}

//=============================================================================
/*
  Create une texture a partir d'un fichier jpg.
  IN: filename: fichier jpg
  IN: mipmap: si true, cree une texture avec mipmaps (not used, always mipmap)
  OUT: tex: texture.
  OUT: width, height: nbre de pixels de la texture
*/
//=============================================================================
E_Int Data::createJpgTexture(const char* filename, GLuint &tex, 
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
  
  struct jpeg_decompress_struct cinfo;
  struct my_error_mgr jerr;
  JSAMPARRAY buffer;		/* Output row buffer */
  int row_stride;		/* physical row width in output buffer */
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = my_error_exit;
  jpeg_create_decompress(&cinfo);

  /* Step 2: specify data source (eg, a file) */
  jpeg_stdio_src(&cinfo, ptrFile);

  /* Step 3: read file parameters with jpeg_read_header() */
  (void) jpeg_read_header(&cinfo, TRUE);
  (void) jpeg_start_decompress(&cinfo);
  row_stride = cinfo.output_width * cinfo.output_components;
  buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

  E_Int nc = cinfo.output_components;
  width = cinfo.output_width;
  height = cinfo.output_height;
  
  bool alpha = false;
  if (nc == 4) alpha = true;

  unsigned char* image = new unsigned char [height * width * nc];

  E_Int j = 0;
  while (cinfo.output_scanline < cinfo.output_height) 
  {
    /* jpeg_read_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could ask for
     * more than one scanline at a time if that's more convenient.
     */
    (void) jpeg_read_scanlines(&cinfo, buffer, 1);
    unsigned char value;
    for (E_Int n = 0; n < nc; n++)
    {
      unsigned char* bufferp = *buffer;
      for (E_Int i = 0; i < width; i++) 
      {
        value = bufferp[nc*i+n];
        image[n+nc*i+(height-j-1)*width*nc] = value;
      }
    }
    j += 1;
  }

  if (tex == 0)
  {
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
    glTexImage2D(GL_TEXTURE_2D, 0, nc, width, height, 0, 
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
    glTexImage2D(GL_TEXTURE_2D, 0, nc, width, height, 0, 
    alpha ? GL_RGBA : GL_RGB, GL_UNSIGNED_BYTE, image);
  }

  fclose(ptrFile);
  return 1;
}
