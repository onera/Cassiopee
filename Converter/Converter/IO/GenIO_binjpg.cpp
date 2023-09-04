/*    
    Copyright 2013-2023 Onera.

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

# include "Images/openjpeg/openjp2/openjpeg.h"
# include "GenIO.h"
# include "Array/Array.h"
# include "Connect/connect.h"
# include <vector>
# include <stdio.h>
# include <string.h>

using namespace K_FLD;
using namespace std;

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

  opj_dparameters_t parameters;               /* Decompression parameters */
  opj_image_t* image = NULL;                  /* Image structure */
  opj_codec_t* l_codec = NULL;                /* Handle to a decompressor */
  opj_stream_t *l_stream = NULL;              /* Stream */

  l_stream = opj_stream_create_default_file_stream(parameters.infile, 1);
  if (!l_stream) 
  {
    fprintf(stderr, "ERROR -> failed to create the stream from the file %s\n",
            parameters.infile);
  }

  switch (parameters.decod_format) 
  {
    case 0: 
    { /* JPEG-2000 codestream */
      /* Get a decoder handle */
      l_codec = opj_create_decompress(OPJ_CODEC_J2K);
      break;
    }

    case 1: 
    { /* JPEG 2000 compressed image data */
      /* Get a decoder handle */
      l_codec = opj_create_decompress(OPJ_CODEC_JP2);
      break;
    }

    case 2: 
    { /* JPEG 2000, JPIP */
      /* Get a decoder handle */
      l_codec = opj_create_decompress(OPJ_CODEC_JPT);
      break;
    }

    default:
    {
      fprintf(stderr, "skipping file..\n");
      //destroy_parameters(&parameters);
      opj_stream_destroy(l_stream);
      return 1;
    }
  }

  //opj_set_info_handler(l_codec, quiet_callback, 00);
  //opj_set_warning_handler(l_codec, quiet_callback, 00);
  //opj_set_error_handler(l_codec, quiet_callback, 00);

  //float t = opj_clock();

  /* Setup the decoder decoding parameters using user parameters */
  //if (!opj_setup_decoder(l_codec, &(parameters.core))) 
  //{
  //  fprintf(stderr, "ERROR -> opj_decompress: failed to setup the decoder\n");
  //  opj_stream_destroy(l_stream);
  //  opj_destroy_codec(l_codec);
  //  return 1;
  //}

  //if (parameters.num_threads >= 1 && !opj_codec_set_threads(l_codec, parameters.num_threads)) 
  //{
  //  fprintf(stderr, "ERROR -> opj_decompress: failed to set number of threads\n");
  //  opj_stream_destroy(l_stream);
  //  opj_destroy_codec(l_codec);
  //  return 1;
  //}

  /* Read the main header of the codestream and if necessary the JP2 boxes */
  if (! opj_read_header(l_stream, l_codec, &image)) 
  {
    fprintf(stderr, "ERROR -> opj_decompress: failed to read the header\n");
    opj_stream_destroy(l_stream);
    opj_destroy_codec(l_codec);
    opj_image_destroy(image);
    return 1;
  }

  // if (parameters.numcomps) 
  // {
  //   if (! opj_set_decoded_components(l_codec,
  //                                     parameters.numcomps,
  //                                     parameters.comps_indices,
  //                                     OPJ_FALSE)) 
  //   {
  //     fprintf(stderr,
  //             "ERROR -> opj_decompress: failed to set the component indices!\n");
  //     opj_destroy_codec(l_codec);
  //     opj_stream_destroy(l_stream);
  //     opj_image_destroy(image);
  //     return 1;
  //   }
  // }

//   if (getenv("USE_OPJ_SET_DECODED_RESOLUTION_FACTOR") != NULL) 
//   {
//     /* For debugging/testing purposes, and also an illustration on how to */
//     /* use the alternative API opj_set_decoded_resolution_factor() instead */
//     /* of setting parameters.cp_reduce */
//     if (! opj_set_decoded_resolution_factor(l_codec, cp_reduce)) 
//     {
//       fprintf(stderr,
//               "ERROR -> opj_decompress: failed to set the resolution factor tile!\n");
//       opj_destroy_codec(l_codec);
//       opj_stream_destroy(l_stream);
//       opj_image_destroy(image);
//       return 1;
//     }
//   }

  if (!parameters.nb_tile_to_decode) 
  {
    if (getenv("SKIP_OPJ_SET_DECODE_AREA") != NULL &&
        parameters.DA_x0 == 0 &&
        parameters.DA_y0 == 0 &&
        parameters.DA_x1 == 0 &&
        parameters.DA_y1 == 0) 
    {
      /* For debugging/testing purposes, */
      /* do nothing if SKIP_OPJ_SET_DECODE_AREA env variable */
      /* is defined and no decoded area has been set */
    }

    /* Optional if you want decode the entire image */
    else if (!opj_set_decode_area(l_codec, image, (OPJ_INT32)parameters.DA_x0,
              (OPJ_INT32)parameters.DA_y0, (OPJ_INT32)parameters.DA_x1,
              (OPJ_INT32)parameters.DA_y1)) 
    {
      fprintf(stderr, "ERROR -> opj_decompress: failed to set the decoded area\n");
      opj_stream_destroy(l_stream);
      opj_destroy_codec(l_codec);
      opj_image_destroy(image);
      return 1;
    }

    /* Get the decoded image */
    if (!(opj_decode(l_codec, l_stream, image) &&
          opj_end_decompress(l_codec, l_stream))) 
    {
      fprintf(stderr, "ERROR -> opj_decompress: failed to decode image!\n");
      opj_destroy_codec(l_codec);
      opj_stream_destroy(l_stream);
      opj_image_destroy(image);
      return 1;
    }
  } 
  else 
  {
    // if (!(parameters.DA_x0 == 0 &&
    //       parameters.DA_y0 == 0 &&
    //       parameters.DA_x1 == 0 &&
    //       parameters.DA_y1 == 0)) 
    // {
    //   if (!(parameters.quiet)) 
    //   {
    //     fprintf(stderr, "WARNING: -d option ignored when used together with -t\n");
    //   }
    // }

    if (!opj_get_decoded_tile(l_codec, l_stream, image, parameters.tile_index)) 
    {
      fprintf(stderr, "ERROR -> opj_decompress: failed to decode tile!\n");
      opj_destroy_codec(l_codec);
      opj_stream_destroy(l_stream);
      opj_image_destroy(image);
      return 1;
    }
    // if (!(parameters.quiet)) 
    // {
    //   fprintf(stdout, "tile %d is decoded!\n\n", parameters.tile_index);
    // }
  }

  /* FIXME? Shouldn't that situation be considered as an error of */
  /* opj_decode() / opj_get_decoded_tile() ? */
  if (image->comps[0].data == NULL) 
  {
    fprintf(stderr, "ERROR -> opj_decompress: no image data!\n");
    opj_destroy_codec(l_codec);
    opj_stream_destroy(l_stream);
    opj_image_destroy(image);
    return 1;
  }

  //tCumulative += opj_clock() - t;
  //numDecompressedImages++;

  /* Close the byte stream */
  opj_stream_destroy(l_stream);

  if (image->color_space != OPJ_CLRSPC_SYCC
      && image->numcomps == 3 && image->comps[0].dx == image->comps[0].dy
      && image->comps[1].dx != 1) 
  {
    image->color_space = OPJ_CLRSPC_SYCC;
  } 
  else if (image->numcomps <= 2) 
  {
    image->color_space = OPJ_CLRSPC_GRAY;
  }

//   if (image->color_space == OPJ_CLRSPC_SYCC) 
//   {
//     color_sycc_to_rgb(image);
//   } 
//   else if ((image->color_space == OPJ_CLRSPC_CMYK) &&
//             (parameters.cod_format != TIF_DFMT)) 
//   {
//     color_cmyk_to_rgb(image);
//   } 
//   else if (image->color_space == OPJ_CLRSPC_EYCC) 
//   {
//     color_esycc_to_rgb(image);
//   }

//   if (image->icc_profile_buf) 
//   {
// #if defined(OPJ_HAVE_LIBLCMS1) || defined(OPJ_HAVE_LIBLCMS2)
//     if (image->icc_profile_len) 
//     {
//       color_apply_icc_profile(image);
//     } 
//     else 
//     {
//       color_cielab_to_rgb(image);
//     }
// #endif
//     free(image->icc_profile_buf);
//     image->icc_profile_buf = NULL;
//     image->icc_profile_len = 0;
//   }

//  /* Force output precision */
//  /* ---------------------- */
//   if (parameters.precision != NULL) 
//   {
//     OPJ_UINT32 compno;
//     for (compno = 0; compno < image->numcomps; ++compno) 
//     {
//       OPJ_UINT32 precno = compno;
//       OPJ_UINT32 prec;

//       if (precno >= parameters.nb_precision) 
//       {
//         precno = parameters.nb_precision - 1U;
//       }

//       prec = parameters.precision[precno].prec;
//       if (prec == 0) 
//       {
//         prec = image->comps[compno].prec;
//       }

//       switch (parameters.precision[precno].mode) 
//       {
//         case OPJ_PREC_MODE_CLIP:
//           clip_component(&(image->comps[compno]), prec);
//           break;
//         case OPJ_PREC_MODE_SCALE:
//           scale_component(&(image->comps[compno]), prec);
//           break;
//         default:
//           break;
//       }

//     }
//   }

//   /* Upsample components */
//   /* ------------------- */
//   if (parameters.upsample) 
//   {
//     image = upsample_image_components(image);
//     if (image == NULL) 
//     {
//       fprintf(stderr,
//               "ERROR -> opj_decompress: failed to upsample image components!\n");
//       opj_destroy_codec(l_codec);
//       return 1;
//     }
//   }

  /* Force RGB output */
  /* ---------------- */
//   if (parameters.force_rgb) 
//   {
//     switch (image->color_space) 
//     {
//       case OPJ_CLRSPC_SRGB:
//         break;
//       case OPJ_CLRSPC_GRAY:
//         image = convert_gray_to_rgb(image);
//         break;
//       default:
//         fprintf(stderr,
//                 "ERROR -> opj_decompress: don't know how to convert image to RGB colorspace!\n");
//         opj_image_destroy(image);
//         image = NULL;
//         break;
//     }
//     if (image == NULL) 
//     {
//       fprintf(stderr, "ERROR -> opj_decompress: failed to convert to RGB image!\n");
//       opj_destroy_codec(l_codec);
//       return 1;
//     }
//   }

  /* create output image */
  /* ------------------- */
//   switch (parameters.cod_format) 
//   {
//     case PXM_DFMT:          /* PNM PGM PPM */
//       if (imagetopnm(image, parameters.outfile, parameters.split_pnm)) 
//       {
//         fprintf(stderr, "[ERROR] Outfile %s not generated\n", parameters.outfile);
//         return 1;
//       } else if (!(parameters.quiet)) 
//       {
//         fprintf(stdout, "[INFO] Generated Outfile %s\n", parameters.outfile);
//       }
//       break;

//       case PGX_DFMT:          /* PGX */
//         if (imagetopgx(image, parameters.outfile)) 
//         {
//           fprintf(stderr, "[ERROR] Outfile %s not generated\n", parameters.outfile);
//           return 1;
//         } 
//         else if (!(parameters.quiet)) 
//         {
//           fprintf(stdout, "[INFO] Generated Outfile %s\n", parameters.outfile);
//         }
//         break;

//         case BMP_DFMT:          /* BMP */
//           if (imagetobmp(image, parameters.outfile)) 
//           {
//             fprintf(stderr, "[ERROR] Outfile %s not generated\n", parameters.outfile);
//             return 1;
//           } 
//           else if (!(parameters.quiet)) 
//           {
//             fprintf(stdout, "[INFO] Generated Outfile %s\n", parameters.outfile);
//           }
//           break;
// #ifdef OPJ_HAVE_LIBTIFF
//         case TIF_DFMT:          /* TIFF */
//             if (imagetotif(image, parameters.outfile)) {
//                 fprintf(stderr, "[ERROR] Outfile %s not generated\n", parameters.outfile);
//                 failed = 1;
//             } else if (!(parameters.quiet)) {
//                 fprintf(stdout, "[INFO] Generated Outfile %s\n", parameters.outfile);
//             }
//             break;
// #endif /* OPJ_HAVE_LIBTIFF */
//         case RAW_DFMT:          /* RAW */
//             if (imagetoraw(image, parameters.outfile)) {
//                 fprintf(stderr, "[ERROR] Error generating raw file. Outfile %s not generated\n",
//                         parameters.outfile);
//                 failed = 1;
//             } else if (!(parameters.quiet)) {
//                 fprintf(stdout, "[INFO] Generated Outfile %s\n", parameters.outfile);
//             }
//             break;

//         case RAWL_DFMT:         /* RAWL */
//             if (imagetorawl(image, parameters.outfile)) {
//                 fprintf(stderr,
//                         "[ERROR] Error generating rawl file. Outfile %s not generated\n",
//                         parameters.outfile);
//                 failed = 1;
//             } else if (!(parameters.quiet)) {
//                 fprintf(stdout, "[INFO] Generated Outfile %s\n", parameters.outfile);
//             }
//             break;

//         case TGA_DFMT:          /* TGA */
//             if (imagetotga(image, parameters.outfile)) {
//                 fprintf(stderr, "[ERROR] Error generating tga file. Outfile %s not generated\n",
//                         parameters.outfile);
//                 failed = 1;
//             } else if (!(parameters.quiet)) {
//                 fprintf(stdout, "[INFO] Generated Outfile %s\n", parameters.outfile);
//             }
//             break;
// #ifdef OPJ_HAVE_LIBPNG
//         case PNG_DFMT:          /* PNG */
//             if (imagetopng(image, parameters.outfile)) {
//                 fprintf(stderr, "[ERROR] Error generating png file. Outfile %s not generated\n",
//                         parameters.outfile);
//                 failed = 1;
//             } else if (!(parameters.quiet)) {
//                 fprintf(stdout, "[INFO] Generated Outfile %s\n", parameters.outfile);
//             }
//             break;
// #endif /* OPJ_HAVE_LIBPNG */
//         /* Can happen if output file is TIFF or PNG
//          * and OPJ_HAVE_LIBTIF or OPJ_HAVE_LIBPNG is undefined
//         */
//         default:
//             fprintf(stderr, "[ERROR] Outfile %s not generated\n", parameters.outfile);
//             failed = 1;
//         }

//         /* free remaining structures */
//         if (l_codec) {
//             opj_destroy_codec(l_codec);
//         }


//         /* free image data structure */
//         opj_image_destroy(image);

//         /* destroy the codestream index */
//         opj_destroy_cstr_index(&cstr_index);
    
  // Stockage du champ
  /*
  varString = new char [128];
  E_Int nil = width;
  E_Int njl = height;
  FldArrayF* f;
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
  */

  opj_stream_destroy(l_stream);
  opj_destroy_codec(l_codec);
  opj_image_destroy(image);
    
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
  vector<E_Int>& eltType,
  vector<char*>& zoneNames)
{
  if (structField.size() == 0) return 0;
  E_Int posR = K_ARRAY::isNamePresent((char*)"r", varString);
  if (posR == -1) posR = K_ARRAY::isNamePresent((char*)"R", varString);
  E_Int posG = K_ARRAY::isNamePresent((char*)"g", varString);
  if (posG == -1) posG = K_ARRAY::isNamePresent((char*)"G", varString);
  E_Int posB = K_ARRAY::isNamePresent((char*)"b", varString);
  if (posB == -1) posB = K_ARRAY::isNamePresent((char*)"B", varString);
  int mode = 0; // only RGB
  
  // recupere la taille de la premiere grille structuree
  int width = 0; int height = 0;
  width = ni[0]; height = nj[0];
  
  //printf("pos %d %d %d - %d %d\n", posR, posG, posB, width, height);
  

  FILE* fp;
  if ((fp = fopen(file, "wb")) == NULL) {
    printf("Warning: cannot open %s.\n", file);
    return 1;
  }
  
  fclose(fp); 

  return 0;
}
