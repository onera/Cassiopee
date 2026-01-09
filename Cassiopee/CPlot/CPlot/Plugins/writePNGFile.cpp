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
#include "Images/png/png.h"
#include "../Data.h"

//============================================================================
// sortie png + depth dans le canal alpha
//============================================================================
void writeDPNGFile(Data* d, char *filename, char* buffer, 
                   E_Int width, E_Int height, E_Int mode)
{
  char* sbuffer = new char [width*height*4];
  for (E_Int j = 0; j < width*height; j++)
  {
    sbuffer[4*j] = buffer[3*j];
    sbuffer[4*j+1] = buffer[3*j+1];
    sbuffer[4*j+2] = buffer[3*j+2];
  }
  for (E_Int j = 0; j < width*height; j++)
  {
    sbuffer[4*j+3] = buffer[3*width*height+j];
  }
  writePNGFile(d, filename, sbuffer, width, height, 1);
  delete [] sbuffer;  
}

//=============================================================================
/*
  Write buffer to PNG file.
  mode=0: RGB, buffer must be RGB
  mode=1: RGB+A, buffer must be RGBA
*/
//=============================================================================
void writePNGFile(Data* d, char* filename, char* buffer, 
                  E_Int width, E_Int height, E_Int mode)
{
  FILE *fp;
  if (buffer == NULL) return;

  png_byte color_type;
  png_byte bit_depth;
  E_Int stride = 0;
  if (mode == 0) 
  { color_type = PNG_COLOR_TYPE_RGB; bit_depth = 8; stride = 3; }
  else { color_type = PNG_COLOR_TYPE_RGBA; bit_depth = 8; stride = 4; }  
  
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep* row_pointers;

  fp = d->fopenw(filename, "wb");
  if (fp == NULL) 
  {
    printf("Warning: cannot open %s.\n", filename);
    return;
  }

  /* initialize stuff */
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_ptr) return;

  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) return;

  if (setjmp(png_jmpbuf(png_ptr))) return;
  png_init_io(png_ptr, fp);

  /* write header */
  if (setjmp(png_jmpbuf(png_ptr))) return;

  png_set_IHDR(png_ptr, info_ptr, width, height,
               bit_depth, color_type, PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

  png_write_info(png_ptr, info_ptr);

  /* write bytes */
  if (setjmp(png_jmpbuf(png_ptr))) return;
   
  row_pointers = new png_bytep [height];
  for (int y = 0; y < height; y++)
  {
    row_pointers[y] = (png_bytep)buffer + (height-1-y)*width*stride;
  }

  png_write_image(png_ptr, row_pointers);

  if (setjmp(png_jmpbuf(png_ptr))) {delete[] row_pointers; return;}
  png_write_end(png_ptr, NULL);

  /* cleanup heap allocation */
  delete [] row_pointers;
  png_destroy_write_struct(&png_ptr, &info_ptr);

  fclose(fp); 

  printf("Wrote file %s (" SF_D_ " x " SF_D_ " pixels).\n", filename, width, height);
}
