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
#include "xcore.h"

#include <stdio.h>
#include <memory.h>
//#include "zoltan/zz/zz_const.h"
//#include "zoltan/zz/zz_rand.h"
//#include "zoltan/params/params_const.h"
//#include "zoltan/all/all_allo_const.h"

/* local function prototypes */
//static void block_part(ZZ *zz, int num_obj, int wtflag, float *wgts, 
//                      float *part_sizes, int *newparts);

/* Block partitioning method.
 * Consider all objects as a linear sequence and do
 * a simple parallel (1d) block partitioning stategy.
 */

// int Zoltan_Block(
//   ZZ *zz,                       /* The Zoltan structure.                     */
//   float *part_sizes,            /* Input:  Array of size zz->LB.Num_Global_Parts
//                                    containing the percentage of work to be
//                                    assigned to each partition.               */
//   int *num_import,              /* Return -1. We use only export lists. */
//   ZOLTAN_ID_PTR *import_global_ids, /* Not used. */
//   ZOLTAN_ID_PTR *import_local_ids,  /* Not used. */
//   int **import_procs,           /* Not used. */
//   int **import_to_part,         /* Not used. */
//   int *num_export,              /* Output: Number of objects to export. */
//   ZOLTAN_ID_PTR *export_global_ids, /* Output: GIDs to export. */
//   ZOLTAN_ID_PTR *export_local_ids,  /* Output: LIDs to export. */
//   int **export_procs,           /* Output: Processsors to export to. */
//   int **export_to_part          /* Output: Partitions to export to. */
// )
// {
//   int ierr = ZOLTAN_OK;
//   int i, count, num_obj;
//   int wtflag = 0;
//   ZOLTAN_ID_PTR global_ids = NULL;
//   ZOLTAN_ID_PTR local_ids = NULL; 
//   int *parts = NULL;
//   int *newparts = NULL;
//   float *wgts = NULL;
//   static char *yo = "Zoltan_Block";

//   ZOLTAN_TRACE_ENTER(zz, yo);

//   /* No import lists computed. */
//   *num_import = -1;
//   *export_global_ids = *export_local_ids = NULL;
//   *export_procs = *export_to_part = NULL;

//   /* Get list of local objects. */
//   if (zz->Obj_Weight_Dim > 1) {
//     ierr = ZOLTAN_FATAL;
//     ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
//                       "OBJ_WEIGHT_DIM > 1 not supported by LB_METHOD BLOCK.");
//     goto End;
//   }
//   wtflag = (zz->Obj_Weight_Dim>0 ? 1 : 0);
//   ierr = Zoltan_Get_Obj_List(zz, &num_obj, &global_ids, &local_ids, wtflag,
//                              &wgts, &parts);

//   /* Compute the new partition numbers. */
//   newparts = (int *) ZOLTAN_MALLOC(num_obj * sizeof(int));
//   if (num_obj && (!newparts)){
//     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
//     ierr = ZOLTAN_MEMERR;
//     goto End;
//   }
//   block_part(zz, num_obj, wtflag, wgts, part_sizes, newparts);

//   /* Check how many partition numbers changed. */
//   count=0;
//   for (i=0; i<num_obj; i++){
//     if (newparts[i] != parts[i])
//       ++count;
//   }
//   (*num_export) = count;

//   /* Allocate export lists. */
//   if ((*num_export) > 0) {
//     if (!Zoltan_Special_Malloc(zz, (void **)export_global_ids, (*num_export),
//                                ZOLTAN_SPECIAL_MALLOC_GID)
//      || !Zoltan_Special_Malloc(zz, (void **)export_local_ids, (*num_export),
//                                ZOLTAN_SPECIAL_MALLOC_LID)
//      || !Zoltan_Special_Malloc(zz, (void **)export_procs, (*num_export),
//                                ZOLTAN_SPECIAL_MALLOC_INT)
//      || !Zoltan_Special_Malloc(zz, (void **)export_to_part, (*num_export),
//                                ZOLTAN_SPECIAL_MALLOC_INT)) {
//       ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
//       ierr = ZOLTAN_MEMERR;
//       goto End;
//     }
//   }

//   /* Loop over objects and fill export lists. */
//   count=0;
//   for (i=0; i<num_obj; i++)
//   {
//     if (newparts[i] != parts[i])
//     {
//       /* export_global_ids[count] = global_ids[i]; */
//       ZOLTAN_SET_GID(zz, &((*export_global_ids)[count*zz->Num_GID]),
//                      &global_ids[i*zz->Num_GID]);
//       if (local_ids)
//         /* export_local_ids[count] = local_ids[i]; */
//         ZOLTAN_SET_LID(zz, &((*export_local_ids)[count*zz->Num_LID]),
//                        &local_ids[i*zz->Num_LID]);
//       /* Set new partition number. */
//       (*export_to_part)[count] = newparts[i];
//       /* Processor is derived from partition number. */
//       (*export_procs)[count] = Zoltan_LB_Part_To_Proc(zz, 
//                      (*export_to_part)[count], &global_ids[i*zz->Num_GID]);
//       ++count;
//     }
//   }

// End:
//   /* Free local memory, but not export lists. */
//   ZOLTAN_FREE(&global_ids);
//   ZOLTAN_FREE(&local_ids);
//   ZOLTAN_FREE(&parts);
//   ZOLTAN_FREE(&newparts);
//   if (wtflag) ZOLTAN_FREE(&wgts);

//   ZOLTAN_TRACE_EXIT(zz, yo);
//   return ierr;
// }


/* Core function to compute partition numbers. */
/* Replace this function with your own algorithm if desired. */
/* Output: newparts contains the new partition numbers. */

// static void block_part(ZZ *zz, int num_obj, int wtflag, float *wgts, 
//             float *part_sizes, int *newparts)
// {
//   int i, part;
//   double wtsum, *scansum;

//   scansum = (double *) ZOLTAN_MALLOC((zz->Num_Proc+1)*sizeof(double));

//   /* Sum up local object weights. */
//   if (wtflag){
//     wtsum = 0.0;
//     for (i=0; i<num_obj; i++) wtsum += wgts[i];
//   }
//   else wtsum = num_obj;

//   /* Cumulative global wtsum */
//   MPI_Allgather(&wtsum, 1, MPI_DOUBLE, &scansum[1], 1, MPI_DOUBLE, 
//                 zz->Communicator);
//   /* scansum = sum of weights on lower processors, excluding self. */
//   scansum[0] = 0.;
//   for (i=1; i<=zz->Num_Proc; i++)
//     scansum[i] += scansum[i-1]; 

//   /* Overwrite part_sizes with cumulative sum (inclusive) part_sizes. */
//   /* A cleaner way is to make a copy, but this works. */
//   for (i=1; i<zz->LB.Num_Global_Parts; i++)
//     part_sizes[i] += part_sizes[i-1]; 

//   /* Loop over objects and assign partition. */
//   part = 0;
//   wtsum = scansum[zz->Proc];
//   for (i=0; i<num_obj; i++){
//     /* wtsum is now sum of all lower-ordered object */
//     /* determine new partition number for this object, 
//        using the "center of gravity" */
//     while (part<zz->LB.Num_Global_Parts-1 && (wtsum+0.5*(wtflag? wgts[i]: 1.0)) 
//            > part_sizes[part]*scansum[zz->Num_Proc])
//       part++;
//     newparts[i] = part;
//     wtsum += (wtflag? wgts[i] : 1.0);
//   }
//   ZOLTAN_FREE(&scansum);
// }

//=============================================================================
// Test 1 de zoltan
//=============================================================================
PyObject* K_XCORE::zoltan1(PyObject* self, PyObject* args)
{
  //PyObject* array;
  //if (!PYPARSETUPLE_(args, O_, &array)) return NULL;

  return NULL;
}
