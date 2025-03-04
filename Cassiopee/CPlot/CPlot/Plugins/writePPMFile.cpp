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
#include "../Data.h"

//=============================================================================
/*
  Write buffer to PPM file.
*/
//=============================================================================
void writePPMFile(Data* d, char *filename, char* buffer, 
                  E_Int width, E_Int height, E_Int mode)
{
  if (buffer == NULL) return;
  E_Int i, j, k, q;
  unsigned char *ibuffer;
  FILE *fp;
  int pixelSize = GL_RGB==GL_RGBA?4:3;
  
  fp = d->fopenw(filename, "wb");
  if ( fp == NULL ) 
  {
    printf("Warning: cannot open %s.\n", filename);
    return;
  }
  
  int elt = 3;
  ibuffer = (unsigned char *)malloc(width*height*3);

  fprintf(fp, "P6\n# CREATOR: CPlot\n" SF_D3_ "\n", width, height, 255);
  q = 0;
  for (i = 0; i < height; i++)
    for (j = 0; j < width; j++)
      for (k = 0; k < elt; k++)
        ibuffer[q++] = (unsigned char)
          *(buffer + (pixelSize*((height-1-i)*width+j)+k));
  fwrite(ibuffer, sizeof(unsigned char), 3*width*height, fp);
  fclose(fp);
  free(ibuffer);
  
  printf("Wrote file %s (" SF_D_ " x " SF_D_ " pixels, " SF_D_ " bytes).\n",
         filename, width, height, 3*width*height);
}
