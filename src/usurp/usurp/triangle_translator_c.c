/* If SINGLE is defined when triangle.o is compiled, it should also be       */
/*   defined here.  If not, it should not be defined here.                   */

/* #define SINGLE */

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

#include <stdio.h>
#include "triangle.h"

#ifndef _STDLIB_H_
extern void *malloc();
extern void free();
#endif /* _STDLIB_H_ */


void report(io, markers, reportpoints, reporttriangles, reportneighbors, 
            reportsegments, reportedges, reportnorms)

struct triangulateio *io;
int markers;
int reportpoints;
int reporttriangles;
int reportneighbors;
int reportsegments;
int reportedges;
int reportnorms;
{
  int i, j;

  if (reportpoints) {
    for (i = 0; i < io->numberofpoints; i++) {
      printf("Point %4d:", i);
      for (j = 0; j < 2; j++) {
        printf("  %.6g", io->pointlist[i * 2 + j]);
      }
      if (io->numberofpointattributes > 0) {
        printf("   attributes");
      }
      for (j = 0; j < io->numberofpointattributes; j++) {
        printf("  %.6g",
               io->pointattributelist[i * io->numberofpointattributes + j]);
      }
      if (markers) {
        printf("   marker %d\n", io->pointmarkerlist[i]);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reporttriangles || reportneighbors) {
    for (i = 0; i < io->numberoftriangles; i++) {
      if (reporttriangles) {
        printf("Triangle %4d points:", i);
        for (j = 0; j < io->numberofcorners; j++) {
          printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
        }
        if (io->numberoftriangleattributes > 0) {
          printf("   attributes");
        }
        for (j = 0; j < io->numberoftriangleattributes; j++) {
          printf("  %.6g", io->triangleattributelist[i *
                                         io->numberoftriangleattributes + j]);
        }
        printf("\n");
      }
      if (reportneighbors) {
        printf("Triangle %4d neighbors:", i);
        for (j = 0; j < 3; j++) {
          printf("  %4d", io->neighborlist[i * 3 + j]);
        }
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reportsegments) {
    for (i = 0; i < io->numberofsegments; i++) {
      printf("Segment %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->segmentlist[i * 2 + j]);
      }
      if (markers) {
        printf("   marker %d\n", io->segmentmarkerlist[i]);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reportedges) {
    for (i = 0; i < io->numberofedges; i++) {
      printf("Edge %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->edgelist[i * 2 + j]);
      }
      if (reportnorms && (io->edgelist[i * 2 + 1] == -1)) {
        for (j = 0; j < 2; j++) {
          printf("  %.6g", io->normlist[i * 2 + j]);
        }
      }
      if (markers) {
        printf("   marker %d\n", io->edgemarkerlist[i]);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }
}

void triangle_translator_c_(int *a, REAL *b)
{
   int ic,iv,i,j,sizea,maxsize,panel;
   struct triangulateio in, out, vorout;


   /* CONVERT THE INTEGER ARRAY */
   i = 0;
   maxsize = a[i++];
   panel = a[i++];
  in.numberofpoints = a[i++];
  in.numberofsegments = a[i++];
  in.numberofholes = a[i++];

  in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
  for (j=0; j < 2*in.numberofsegments; j++) {
      in.segmentlist[j] = a[i++];
   }

  /* CONVERT THE REAL ARRAY */
  in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
  i = 0;
  for (j=0; j < 2*in.numberofpoints; j++) {
      in.pointlist[j] = b[i++];
   }


  /* other stuff we will not need */
  in.numberofregions = 0;
  in.regionlist = (REAL *) NULL;
  in.numberofpointattributes = 0;
  in.pointattributelist = (REAL *) NULL;
  in.pointmarkerlist = (int *) NULL;
  in.segmentmarkerlist = (int *) NULL;

  in.numberoftriangles = 0;
  in.numberofcorners = 3;
  in.numberoftriangleattributes = 0;
  in.trianglelist = (int *) NULL;
  in.triangleattributelist = (REAL *) NULL;
  in.trianglearealist = (REAL *) NULL;
  in.holelist = (REAL *) NULL;

  /* printf("Input point set:\n\n");
  report(&in, 0, 1, 0, 0, 0, 0, 0); */

  /* Make necessary initializations so that Triangle can return a */
  /*   triangulation in `out'                                     */

  out.trianglelist = (int *) NULL;

  /* Triangulate the points.  Switches are chosen to  */
  /*   (p) read and write a PSLG                      */
  /*   (z) number everything from zero                */
  /*   (N) do not output nodes (same as input)        */
  /*   (B) no boundary markers in output              */
  /*   (Q) quiet                                      */
  /*   (P) do not output poly file                    */

  triangulate("pzBPNQ", &in, &out, &vorout);

/*   if (panel == 16417) {
     printf("Triangulation:\n\n");
     report(&out, 0, 0, 1, 0, 0, 0, 0); 
   } */


   /* COUNT THE INTEGERS FROM THE RESULT */
   /* include total number of integers in data stream */
   /* include total number of triangles */
   /* include number of corners (expect 3) */
   /* include actual corner data (number_corners * number_triangles) */

   sizea = 3 + out.numberoftriangles * out.numberofcorners;
   if (sizea > maxsize) {
      a[0] = -1*sizea;

      /* deallocate data */
      free(in.pointlist);
      free(in.segmentlist);
      free(out.trianglelist);
      return;
   }

   /* STORE THE INTEGERS FROM THE RESULT */
    i = 0;
    a[i++] = sizea;
    a[i++] = out.numberoftriangles;
    a[i++] = out.numberofcorners;
    for (ic = 0; ic < out.numberoftriangles; ic++) {
      for (iv = 0; iv < out.numberofcorners; iv++) {
        a[i++] = out.trianglelist[ic * out.numberofcorners + iv];
      }
    }


   /* COUNT THE FLOATS FROM THE RESULT */
   /* not needed for now since output nodes same as input nodes */

   /* STORE THE FLOATS FROM THE RESULT */
   /* not needed for now since output nodes same as input nodes */


  /* Free all allocated arrays, including those allocated by Triangle. */

  free(in.pointlist);
  free(in.segmentlist);
  free(out.trianglelist);

  return;

}
