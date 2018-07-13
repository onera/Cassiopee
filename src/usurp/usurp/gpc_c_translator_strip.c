#include <stdlib.h>
#include "gpc.h"

void gpc_c_translator_strip_(int *a, double *b)
{
   int ic,iv,i,sizea,sizeb,maxsize;
   gpc_polygon subj;
   gpc_tristrip tristrip;


   /* CONVERT THE INTEGER ARRAY */
   i = 0;
   maxsize = a[i++];

   subj.num_contours = a[i++];

   subj.contour = malloc(subj.num_contours * sizeof(gpc_vertex_list));
   subj.hole    = malloc(subj.num_contours * sizeof(int));

   for (ic=0; ic < subj.num_contours; ic++) {
       subj.contour[ic].num_vertices = a[i++];
       subj.hole[ic] = a[i++];
   }

   /* CONVERT THE FLOAT ARRAY */
   i = 0;
   for (ic=0; ic < subj.num_contours; ic++) {

       subj.contour[ic].vertex = malloc(subj.contour[ic].num_vertices * sizeof(gpc_vertex));

       for (iv = 0; iv < subj.contour[ic].num_vertices; iv++) { 
          subj.contour[ic].vertex[iv].x = b[i++];
          subj.contour[ic].vertex[iv].y = b[i++];
       }
   } 


   /* CALL THE GPC LIBRARY */
   gpc_polygon_to_tristrip(&subj, &tristrip);


   /* COUNT THE INTEGERS FROM THE RESULT */
    sizea = 2 + 2 * tristrip.num_strips ;
    if (sizea > maxsize) {
       a[0] = -1*sizea;

       gpc_free_polygon(&subj);
       gpc_free_tristrip(&tristrip);
       return;
    }


   /* STORE THE INTEGERS FROM THE RESULT */
    i = 0;
    a[i++] = sizea;
    a[i++] = tristrip.num_strips;
    for (ic = 0; ic < tristrip.num_strips; ic++) {
       a[i++] = tristrip.strip[ic].num_vertices;
       a[i++] = 0;
    }
    
   /* COUNT THE FLOATS FROM THE RESULT */
    sizeb = 0;
    for (ic = 0; ic < tristrip.num_strips; ic++) {
      sizeb = sizeb + 2 * tristrip.strip[ic].num_vertices;
    }
    if (sizeb > maxsize)
    {
       a[0] = -1*sizeb;

       gpc_free_polygon(&subj);
       gpc_free_tristrip(&tristrip);
       return;
    }

   /* STORE THE FLOATS FROM THE RESULT */
    i = 0;
    for (ic = 0; ic < tristrip.num_strips; ic++) {
      for (iv = 0; iv < tristrip.strip[ic].num_vertices; iv++) {
          b[i++] = tristrip.strip[ic].vertex[iv].x;
          b[i++] = tristrip.strip[ic].vertex[iv].y;
      }
    }

   gpc_free_polygon(&subj);
   gpc_free_tristrip(&tristrip);
   return;
}
