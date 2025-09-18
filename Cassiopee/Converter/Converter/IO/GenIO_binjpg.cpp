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

// Binary JPG (JPEG) file support

# include "GenIO.h"
# include "Array/Array.h"
# include "String/kstring.h"
# include "Connect/connect.h"
# include <vector>
# include <stdio.h>
# include <string.h>

#include "Images/libjpeg/jpeglib.h"

using namespace K_FLD;
using namespace std;

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
   jpgread 
   Un seul domaine peut-etre stocke dans ce type de fichier.
*/
//=============================================================================
E_Int K_IO::GenIO::jpgread( 
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
    printf("Warning: jpgread: cannot open file %s.\n", file);
    return 1;
  }

  /* This struct contains the JPEG decompression parameters and pointers to
   * working space (which is allocated as needed by the JPEG library).
   */
  struct jpeg_decompress_struct cinfo;
  /* We use our private extension JPEG error handler.
   * Note that this struct must live as long as the main JPEG parameter
   * struct, to avoid dangling-pointer problems.
   */
  struct my_error_mgr jerr;
  /* More stuff */
  JSAMPARRAY buffer;		/* Output row buffer */
  int row_stride;		/* physical row width in output buffer */

  /* In this example we want to open the input file before doing anything else,
   * so that the setjmp() error recovery below can assume the file is open.
   * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
   * requires it in order to read binary files.
   */
  
  /* Step 1: allocate and initialize JPEG decompression object */

  /* We set up the normal JPEG error routines, then override error_exit. */
  
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = my_error_exit;
  /* Establish the setjmp return context for my_error_exit to use. */
  //if (setjmp(jerr.setjmp_buffer)) {
    /* If we get here, the JPEG code has signaled an error.
     * We need to clean up the JPEG object, close the input file, and return.
     */
  //  jpeg_destroy_decompress(&cinfo);
  //  fclose(ptrFile);
  //  return 0;
  //}

  /* Now we can initialize the JPEG decompression object. */
  jpeg_create_decompress(&cinfo);

  /* Step 2: specify data source (eg, a file) */
  jpeg_stdio_src(&cinfo, ptrFile);

  /* Step 3: read file parameters with jpeg_read_header() */
  (void) jpeg_read_header(&cinfo, TRUE);

  /* We can ignore the return value from jpeg_read_header since
   *   (a) suspension is not possible with the stdio data source, and
   *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
   * See libjpeg.txt for more info.
   */

  /* Step 4: set parameters for decompression */

  /* In this example, we don't need to change any of the defaults set by
   * jpeg_read_header(), so we do nothing here.
   */

  /* Step 5: Start decompressor */
  (void) jpeg_start_decompress(&cinfo);
  
  /* We can ignore the return value since suspension is not possible
   * with the stdio data source.
   */
  /* We may need to do some setup of our own at this point before reading
   * the data.  After jpeg_start_decompress() we have the correct scaled
   * output image dimensions available, as well as the output colormap
   * if we asked for color quantization.
   * In this example, we need to make an output work buffer of the right size.
   */ 
  /* JSAMPLEs per row in output buffer */
  row_stride = cinfo.output_width * cinfo.output_components;
  /* Make a one-row-high sample array that will go away when done with image */
  buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

  // allocate output Fld
  FldArrayF* f;
  E_Int nc = cinfo.output_components;
  E_Int nil = cinfo.output_width;
  E_Int njl = cinfo.output_height;
  //printf("size=" SF_D2_ ", components=" SF_D_ "\n", nil, njl, nc);
  
  varString = new char [128];

  if (nc == 1)
  {
    strcpy(varString, "x,y,z,r");
    f = new FldArrayF(nil*njl, 4);
  }
  else if (nc == 2)
  {
    strcpy(varString, "x,y,z,r,g");
    f = new FldArrayF(nil*njl, 5);
  }
  else if (nc == 3)
  {
    strcpy(varString, "x,y,z,r,g,b");
    f = new FldArrayF(nil*njl, 6);
  }
  else if (nc == 4)
  {
    strcpy(varString, "x,y,z,r,g,b,a");
    f = new FldArrayF(nil*njl, 7);
  }
  else
  {
    printf("readjpg: the number of components is not possible.");
    exit(0);
  }

  f->setAllValuesAtNull();
  E_Float* fx = f->begin(1);
  E_Float* fy = f->begin(2);
  E_Float* fz = f->begin(3);

  // Cree les noms de zones
  char* zoneName = new char [128];
  sprintf(zoneName, "Image0");
  zoneNames.push_back(zoneName);

  for (E_Int j = 0; j < njl; j++)
    for (E_Int i = 0; i < nil; i++)
    {
      fx[i+j*nil] = i*1.;
      fy[i+j*nil] = (njl-j)*1.;
      fz[i+j*nil] = 0.;
    }

  structField.push_back(f);
  ni.push_back(nil); nj.push_back(njl); nk.push_back(1);

  /* Step 6: while (scan lines remain to be read) */
  /*           jpeg_read_scanlines(...); */

  /* Here we use the library's state variable cinfo.output_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   */
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
      E_Float* fr = f->begin(n+4);
      unsigned char* bufferp = *buffer;
      for (E_Int i = 0; i < nil; i++) 
      {
        value = bufferp[nc*i+n];
        fr[i+j*nil] = (E_Float)value;
      }
    }
    j += 1;
  }

  /* Step 7: Finish decompression */
  (void) jpeg_finish_decompress(&cinfo);
  /* We can ignore the return value since suspension is not possible
   * with the stdio data source.
   */

  /* Step 8: Release JPEG decompression object */

  /* This is an important step since it will release a good deal of memory. */
  jpeg_destroy_decompress(&cinfo);
    
  fclose(ptrFile);
  return 0;
}

//=============================================================================
// Ecrit seulement le premier domaine structure
// Recupere seulement les champs RGB
// Les champs doivent etre en entiers entre 0 et 255
//=============================================================================
E_Int K_IO::GenIO::jpgwrite(
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
  

  //printf("pos " SF_D3_ " - " SF_D2_ "\n", posR, posG, posB, width, height);
  E_Int nc = 3; // only for RGB for now

  FILE* ptrFile;
  if ((ptrFile = fopen(file, "wb")) == NULL) {
    printf("Warning: cannot open %s.\n", file);
    return 1;
  }
  
  /* passe le Fld dans image_buffer */
  JSAMPLE* image_buffer;	/* Points to large array of R,G,B-order data */
  int image_width = ni[0];	/* Number of rows in image */
  int image_height = nj[0];		/* Number of columns in image */
  image_buffer = new JSAMPLE [image_width*nc*image_height];
  E_Float* r = structField[0]->begin(posR+1);
  E_Float* g = structField[0]->begin(posG+1);
  E_Float* b = structField[0]->begin(posB+1);
  
  JSAMPLE* p = image_buffer;
  for (E_Int j = 0; j < image_height; j++)
  {
    for (E_Int i = 0; i < image_width; i++)
    {
      *p = (unsigned char)r[i+j*image_width]; p++;
      *p = (unsigned char)g[i+j*image_width]; p++;
      *p = (unsigned char)b[i+j*image_width]; p++;
    }
  }

  struct jpeg_compress_struct cinfo;
  
  struct jpeg_error_mgr jerr;
  JSAMPROW row_pointer[1];	/* pointer to JSAMPLE row[s] */
  int row_stride;		/* physical row width in image buffer */

  cinfo.err = jpeg_std_error(&jerr);
  /* Now we can initialize the JPEG compression object. */
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, ptrFile);

  cinfo.image_width = image_width; 	/* image width and height, in pixels */
  cinfo.image_height = image_height;
  cinfo.input_components = 3;		/* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
  /* Now use the library's routine to set default compression parameters.
   * (You must set at least cinfo.in_color_space before calling this,
   * since the defaults depend on the source color space.)
   */
  jpeg_set_defaults(&cinfo);
  /* Now you can set any non-default parameters you wish to.
   * Here we just illustrate the use of quality (quantization table) scaling:
   */
  int quality = 100; // a ajuster
  jpeg_set_quality(&cinfo, quality, TRUE /* limit to baseline-JPEG values */);

  /* Step 4: Start compressor */

  /* TRUE ensures that we will write a complete interchange-JPEG file.
   * Pass TRUE unless you are very sure of what you're doing.
   */
  jpeg_start_compress(&cinfo, TRUE);

  /* Step 5: while (scan lines remain to be written) */
  /*           jpeg_write_scanlines(...); */

  /* Here we use the library's state variable cinfo.next_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   * To keep things simple, we pass one scanline per call; you can pass
   * more if you wish, though.
   */
  row_stride = image_width * 3;	/* JSAMPLEs per row in image_buffer */

  while (cinfo.next_scanline < cinfo.image_height) {
    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    row_pointer[0] = & image_buffer[cinfo.next_scanline * row_stride];
    (void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }

  /* Step 6: Finish compression */

  jpeg_finish_compress(&cinfo);
  /* After finish_compress, we can close the output file. */
  fclose(ptrFile);

  /* Step 7: release JPEG compression object */

  /* This is an important step since it will release a good deal of memory. */
  jpeg_destroy_compress(&cinfo);

  return 0;
}
