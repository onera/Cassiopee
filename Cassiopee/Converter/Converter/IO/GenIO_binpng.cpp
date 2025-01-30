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

// Binary PNG (Portable Network Graphics) file support

# include "Images/png/png.h"
# include "GenIO.h"
# include "Array/Array.h"
# include "String/kstring.h"
# include "Connect/connect.h"
# include <vector>
# include <stdio.h>
# include <string.h>

using namespace K_FLD;
using namespace std;

//=============================================================================
/* 
   pngread 
   Un seul domaine peut-etre stocke dans ce type de fichier.
   Limite a un bit_depth de 8.
*/
//=============================================================================
E_Int K_IO::GenIO::pngread( 
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, 
  vector<char*>& zoneNames)
{
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "rb");

  if (ptrFile == NULL)
  {
    printf("Warning: pngread: cannot open file %s.\n", file);
    return 1;
  }

  int width, height;
  png_structp png_ptr;
  png_infop info_ptr, end_info_ptr;
  //int number_of_passes;
  png_bytep* row_pointers;
  png_byte color_type;
  png_byte bit_depth;

  // Header
  png_byte header[8];    // 8 is the maximum size that can be checked
  fread(header, 1, 8, ptrFile);
  if (png_sig_cmp(header, 0, 8))
  { printf("File %s is an invalid png file.\n", file); fclose(ptrFile); return 1; }
  
  /* initialize stuff */
  png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (!png_ptr)
  { printf("Can not create png structures.\n"); fclose(ptrFile); return 1; }
  
  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr)
  { printf("Can not create png structures.\n"); fclose(ptrFile); return 1; }
               
  end_info_ptr = png_create_info_struct(png_ptr);
  if (!end_info_ptr)
  { printf("Can not create png structures.\n"); fclose(ptrFile); return 1; }
 
  if (setjmp(png_jmpbuf(png_ptr)))
  { printf("Can not create png structures.\n"); fclose(ptrFile); return 1; }
              
  png_init_io(png_ptr, ptrFile);
  png_set_sig_bytes(png_ptr, 8);

  png_read_info(png_ptr, info_ptr);

  width = png_get_image_width(png_ptr, info_ptr);
  height = png_get_image_height(png_ptr, info_ptr);
  color_type = png_get_color_type(png_ptr, info_ptr);
  bit_depth = png_get_bit_depth(png_ptr, info_ptr);
  //printf("color type " SF_D_ "\n", color_type);
  //printf("bit depth " SF_D_ "\n", bit_depth);
  
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

  //number_of_passes = png_set_interlace_handling(png_ptr);
  png_read_update_info(png_ptr, info_ptr);

  varString = new char [128];

  /* read file */
  if (setjmp(png_jmpbuf(png_ptr)))
  { printf("Error during read of png image.\n"); fclose(ptrFile); return 1; }

  row_pointers = (png_bytep*)malloc(sizeof(png_bytep)*height);
  for (int y = 0; y < height; y++)
    row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png_ptr, info_ptr));

  png_read_image(png_ptr, row_pointers);

  int components;
  int size = png_get_rowbytes(png_ptr, info_ptr);
  components = (size/width);
  //printf("components = " SF_D_ "\n", components);

  // Stockage du champ
  E_Int nil = width;
  E_Int njl = height;
  FldArrayF* f = NULL;
  if (components == 1) // greyscale
  {
    strcpy(varString, "x,y,z,r");
    f = new FldArrayF(nil*njl, 4);
  }
  else if (components == 3) // RGB
  {
    strcpy(varString, "x,y,z,r,g,b");
    f = new FldArrayF(nil*njl, 6);
  }
  else if (components == 2) // greyscale + alpha
  {
    strcpy(varString, "x,y,z,r,a");
    f = new FldArrayF(nil*njl, 5);
  }
  else if (components == 4) // RGB + alpha
  {
    strcpy(varString, "x,y,z,r,g,b,a");
    f = new FldArrayF(nil*njl, 7);
  }
  else
  {
    printf("Warning: bin_png: unkown number of components.\n");
  }
 
  f->setAllValuesAtNull();
  E_Float* fx = f->begin(1);
  E_Float* fy = f->begin(2);
  E_Float* fz = f->begin(3);

  structField.push_back(f);
  ni.push_back(nil); nj.push_back(njl); nk.push_back(1);
  
  // Cree les noms de zones
  char* zoneName = new char [128];
  sprintf(zoneName, "Zone0");
  zoneNames.push_back(zoneName);

  for (E_Int j = 0; j < njl; j++)
    for (E_Int i = 0; i < nil; i++)
    {
      fx[i+j*nil] = i*1.;
      fy[i+j*nil] = (njl-j)*1.;
      fz[i+j*nil] = 0.;
    }

  if (components == 1) // greyscale
  {
    E_Float* r = f->begin(4);
    for (E_Int j = 0; j < njl; j++)
    {
      png_byte* p = row_pointers[j];
      for (E_Int i = 0; i < nil; i++)
      {
        r[i+j*nil] = *p; p++;
      }
    }
  }
  else if (components == 3)
  {
    E_Float* r = f->begin(4);
    E_Float* g = f->begin(5);
    E_Float* b = f->begin(6);
    for (E_Int j = 0; j < njl; j++)
    {
      png_byte* p = row_pointers[j];
      for (E_Int i = 0; i < nil; i++)
      {
        r[i+j*nil] = *p; p++;
        g[i+j*nil] = *p; p++;
        b[i+j*nil] = *p; p++;
      }
    }
  }
  else if (components == 2)
  {
    E_Float* r = f->begin(4);
    E_Float* alpha = f->begin(5);
    for (E_Int j = 0; j < njl; j++)
    {
      png_byte* p = row_pointers[j];
      for (E_Int i = 0; i < nil; i++)
      {
        r[i+j*nil] = *p; p++;
        alpha[i+j*nil] = *p; p++;
      }
    }
  }
  else if (components == 4)
  {
    E_Float* r = f->begin(4);
    E_Float* g = f->begin(5);
    E_Float* b = f->begin(6);
    E_Float* alpha = f->begin(7);
    for (E_Int j = 0; j < njl; j++)
    {
      png_byte* p = row_pointers[j];
      for (E_Int i = 0; i < nil; i++)
      {
        r[i+j*nil] = *p; p++;
        g[i+j*nil] = *p; p++;
        b[i+j*nil] = *p; p++;
        alpha[i+j*nil] = *p; p++;
      }
    }
  }

  for (int y = 0; y < height; y++) free(row_pointers[y]);
  free(row_pointers);
  
  png_read_end(png_ptr, end_info_ptr);
  png_destroy_read_struct(&png_ptr, &info_ptr, &end_info_ptr);

  fclose(ptrFile);
  return 0;
}

//=============================================================================
// Ecrit seulement le premier domaine structure
// Recupere seulement les champs RGB
// Les champs doivent etre en entiers entre 0 et 255
//=============================================================================
E_Int K_IO::GenIO::pngwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector< vector<E_Int> >& eltTypes,
  vector<char*>& zoneNames)
{
  if (structField.size() == 0) return 0;
  E_Int posR = K_ARRAY::isNamePresent((char*)"r", varString);
  if (posR == -1) posR = K_ARRAY::isNamePresent((char*)"R", varString);
  E_Int posG = K_ARRAY::isNamePresent((char*)"g", varString);
  if (posG == -1) posG = K_ARRAY::isNamePresent((char*)"G", varString);
  E_Int posB = K_ARRAY::isNamePresent((char*)"b", varString);
  if (posB == -1) posB = K_ARRAY::isNamePresent((char*)"B", varString);
  E_Int posA = K_ARRAY::isNamePresent((char*)"a", varString);
  if (posA == -1) posA = K_ARRAY::isNamePresent((char*)"A", varString);
  
  int mode = 0; // only RGB
  
  // recupere la taille de la premiere grille structuree
  int width = 0; int height = 0;
  width = ni[0]; height = nj[0];
  
  //printf("pos " SF_D3_ " - " SF_D2_ "\n", posR, posG, posB, width, height);
  
  E_Int nc = 3;
  if (posA >= 0) { nc = 4; mode = 1; }

  // cree le buffer
  png_byte* buffer = new png_byte [nc*width*height];
  for (E_Int i = 0; i < nc*width*height; i++) buffer[i] = 0;
  if (posR >= 0)
  {
    E_Float* r = structField[0]->begin(posR+1);
    for (E_Int i = 0; i < width*height; i++) buffer[nc*i] = (png_byte)r[i];
  }
  if (posG >= 0)
  {
    E_Float* g = structField[0]->begin(posG+1);
    for (E_Int i = 0; i < width*height; i++) buffer[nc*i+1] = (png_byte)g[i];
  }
  if (posB >= 0)
  {
    E_Float* b = structField[0]->begin(posB+1);
    for (E_Int i = 0; i < width*height; i++) buffer[nc*i+2] = (png_byte)b[i];
  }
  if (posA >= 0)
  {
    E_Float* a = structField[0]->begin(posA+1);
    for (E_Int i = 0; i < width*height; i++) buffer[nc*i+3] = (png_byte)a[i];
  }
  
  png_byte color_type;
  png_byte bit_depth;
  int stride = 0;
  if (mode == 0)
  { color_type = PNG_COLOR_TYPE_RGB; bit_depth = 8; stride = 3; }
  else { color_type = PNG_COLOR_TYPE_RGBA; bit_depth = 8; stride = 4; }  
  
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep* row_pointers;

  FILE* fp;
  if ((fp = fopen(file, "wb")) == NULL) {
    printf("Warning: cannot open %s.\n", file);
    return 1;
  }

  /* initialize stuff */
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_ptr) return 1;

  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) return 1;

  if (setjmp(png_jmpbuf(png_ptr))) return 1;
  png_init_io(png_ptr, fp);

  /* write header */
  if (setjmp(png_jmpbuf(png_ptr))) return 1;

  png_set_IHDR(png_ptr, info_ptr, width, height,
               bit_depth, color_type, PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

  png_write_info(png_ptr, info_ptr);

  /* write bytes */
  if (setjmp(png_jmpbuf(png_ptr))) return 1;
  
  row_pointers = new png_bytep [height];
  for (int y = 0; y < height; y++)
  {
    //row_pointers[y] = (png_bytep)buffer + (height-1-y)*width*stride;
    row_pointers[y] = (png_bytep)buffer + y*width*stride;
  }

  png_write_image(png_ptr, row_pointers);

  if (setjmp(png_jmpbuf(png_ptr))) {delete[] row_pointers; return 1;}

  png_write_end(png_ptr, NULL);
  png_destroy_write_struct(&png_ptr, &info_ptr);

  /* cleanup allocations */
  delete [] row_pointers;
  delete [] buffer;
  
  fclose(fp); 

  return 0;
}
