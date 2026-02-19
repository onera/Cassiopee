#include <stdlib.h>
#include "gpc.h"

extern "C"
{
  void gpc_c_translator_clip_(int *a, double *b);
}

void gpc_c_translator_clip_(int *a, double *b)
{
   int ic,iv,i,sizea,sizeb,maxsize;
   gpc_polygon subj, clip, result;


   /* CONVERT THE INTEGER ARRAY */
   i = 0;
   maxsize = a[i++];

   subj.num_contours = a[i++];

   subj.contour = (gpc_vertex_list *)
     malloc(subj.num_contours * sizeof(gpc_vertex_list));
   subj.hole    = (int*)malloc(subj.num_contours * sizeof(int));

   for (ic=0; ic < subj.num_contours; ic++) 
   {
       subj.contour[ic].num_vertices = a[i++];
       subj.hole[ic] = a[i++];
   }

   clip.num_contours = a[i++];

   clip.contour = (gpc_vertex_list *)
     malloc(clip.num_contours * sizeof(gpc_vertex_list));
   clip.hole = (int*)malloc(clip.num_contours * sizeof(int));

   for (ic=0; ic < clip.num_contours; ic++) {
       clip.contour[ic].num_vertices = a[i++];
       clip.hole[ic] = a[i++];
   }

   /* CONVERT THE FLOAT ARRAY */
   i = 0;
   for (ic=0; ic < subj.num_contours; ic++) {

     subj.contour[ic].vertex = (gpc_vertex*)
       malloc(subj.contour[ic].num_vertices * sizeof(gpc_vertex));

       for (iv = 0; iv < subj.contour[ic].num_vertices; iv++) { 
          subj.contour[ic].vertex[iv].x = b[i++];
          subj.contour[ic].vertex[iv].y = b[i++];
       }
   } 
   for (ic=0; ic < clip.num_contours; ic++) {

       clip.contour[ic].vertex = 
         (gpc_vertex*)malloc(clip.contour[ic].num_vertices * sizeof(gpc_vertex));

       for (iv = 0; iv < clip.contour[ic].num_vertices; iv++) { 
          clip.contour[ic].vertex[iv].x = b[i++];
          clip.contour[ic].vertex[iv].y = b[i++];
       }
   }


   /* SANITY CHECKS */
/* if (clip.num_contours != 1) {
      printf("ERROR! clip.num_contours = %d\n",clip.num_contours);
      a[0] = -a[0];
      return;
   }
   if (clip.contour[0].num_vertices != 4) {
      printf("ERROR! clip.contour[1].num_vertices = %d\n",clip.contour[1].num_vertices);
      printf("ERROR! clip.num_contours = %d\n",clip.num_contours);
      i = 0;

      printf("a[%d] = %d (max size)\n",i,a[i]);
      i++;
      printf("maxsize           = %d\n",maxsize);

      printf("a[%d] = %d (subj.num_contours)\n",i,a[i]);
      i++;
      printf("subj.num_contours = %d\n",subj.num_contours);

      for (ic = 0; ic < subj.num_contours; ic++) {
          printf("a[%d] = %d (subj.num_verts)\n",i,a[i]);
          i++;
          printf("subj.contour[%d].num_vertices = %d\n",ic,subj.contour[ic].num_vertices);
          printf("a[%d] = %d (subj.hole)\n",i,a[i]);
          i++;
      }

      printf("a[%d] = %d (clip.num_contours)\n",i,a[i]);
      i++;
      printf("clip.num_contours = %d\n",clip.num_contours);

      for (ic = 0; ic < clip.num_contours; ic++) {
          printf("a[%d] = %d (clip.num_verts)\n",i,a[i]);
          i++;
          printf("clip.contour[%d].num_vertices = %d\n",ic,clip.contour[ic].num_vertices);
          printf("a[%d] = %d (clip.hole)\n",i,a[i]);
          i++;
      }

      a[0] = -a[0];
      return;
   } */

   /* CALL THE GPC LIBRARY */
   gpc_polygon_clip(GPC_DIFF, &subj, &clip, &result); 
   if (result.num_contours < 0){
      fprintf(stderr," ERROR! gpc_polygon_clip returned result.num_contours < 0\n");
      a[0] = result.num_contours;
      gpc_free_polygon(&subj);
      gpc_free_polygon(&clip);
      gpc_free_polygon(&result);
      return;
   }

   /* COUNT THE INTEGERS FROM THE RESULT */
   sizea = 2 + 2 * result.num_contours ;
   if (sizea > maxsize) {
      a[0] = -1*sizea;

      gpc_free_polygon(&subj);
      gpc_free_polygon(&clip);
      gpc_free_polygon(&result);
      return;
   }

   /* CHECK THE HOLE VALUES */
   for (ic = 0; ic < result.num_contours; ic++) {
      if ((result.hole[ic] != 0) &&  (result.hole[ic] != 1)) {
        a[0] = -2;
        gpc_free_polygon(&subj);
        gpc_free_polygon(&clip);
        gpc_free_polygon(&result);
        return;
      }
   }
      

   /* STORE THE INTEGERS FROM THE RESULT */
   i = 0;
   a[i++] = sizea;
   a[i++] = result.num_contours;
   for (ic = 0; ic < result.num_contours; ic++) {
      a[i++] = result.contour[ic].num_vertices;
      a[i++] = result.hole[ic];
   }
 

   /* COUNT THE FLOATS FROM THE RESULT */
   sizeb = 0;
   for (ic = 0; ic < result.num_contours; ic++) {
      sizeb = sizeb + 2 * result.contour[ic].num_vertices;
   }
   if (sizeb > maxsize)
   {
      a[0] = -1*sizeb;

      gpc_free_polygon(&subj);
      gpc_free_polygon(&clip);
      gpc_free_polygon(&result);
      return;
   }


   /* STORE THE FLOATS FROM THE RESULT */
   i = 0;
   for (ic = 0; ic < result.num_contours; ic++) {
      for (iv = 0; iv < result.contour[ic].num_vertices; iv++) {
         b[i++] = result.contour[ic].vertex[iv].x;
         b[i++] = result.contour[ic].vertex[iv].y;
     }
   }

   gpc_free_polygon(&subj);
   gpc_free_polygon(&clip);
   gpc_free_polygon(&result);
   return;
}
