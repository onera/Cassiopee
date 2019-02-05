
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_geom_elem.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Maximum number of sections depending of section type
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_BLOCK_ID_BLOCK_STD    = 0,     
  PDM_BLOCK_ID_BLOCK_POLY2D = 1000000,
  PDM_BLOCK_ID_BLOCK_POLY3D = 2000000

} PDM_section_id_section_t;

/*============================================================================
 * Global variable
 *============================================================================*/

/**
 * \brief Storage of mesh handles
 *
 */

static PDM_Handles_t *mesh_handles = NULL; 

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Quick sort
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 * 
 */

static void 
_quickSort_int
( 
 int a[], 
 int l, 
 int r 
)
{
  if (l < r) {
    int j = r+1;
    int t; 
    int pivot = a[l];
    int i = l; 

    while(1) {
      do ++i; while (a[i] <= pivot && i < r);
      do --j; while (a[j] > pivot);
      if (i >= j) break;

      t    = a[i]; 
      a[i] = a[j]; 
      a[j] = t;

    }
    t    = a[l]; 
    a[l] = a[j]; 
    a[j] = t;

    _quickSort_int(a, l  , j-1);
    _quickSort_int(a, j+1,   r);
  }
}

/**
 * \def _find_pairs
 * Search common faces in a distribution 
 *
 */
static int
_find_pairs
(
const int *IdxFace,
const int *data,
const int  nFac,
      int  iAbsFace,  /* A passer ren reference ou a return */
      PDM_g_num_t *dFaceVtx,
      int *dFaceVtxIdx,
      PDM_g_num_t *dFaceCell
)
{
  if(0 == 1){
    printf("_find_pairs : nFac : %d \n", nFac);
  }

  /*
   * Special Case -> To removed when Boundary is OK
   */
  // if(nFac == 1){
  //   // printf("Something strange append for in find_pairs - Check mesh ooo \n");
    
  //   /* Stokage en element externe */
  //   int iFacIdx = dFaceVtxIdx[iAbsFace];
    
  //   int curFac = IdxFace[0];
  //   int nVtx1  = data[curFac+2];
    
    
  //   dFaceCell[2*iAbsFace  ] = data[curFac];
  //   dFaceCell[2*iAbsFace+1] = 0;
    
  //   // printf("----------: %d - %d - %d\n", curFac, nVtx1, iFacIdx);
 
    
  //   for(int iVtx=0; iVtx < nVtx1; iVtx++){
  //     dFaceVtx[iFacIdx+iVtx] = data[curFac+3+iVtx];
  //   }
    
  //   /* Update index */
  //   dFaceVtxIdx[iAbsFace+1] = dFaceVtxIdx[iAbsFace] + nVtx1;
  //   iAbsFace++;
    
  //   return iAbsFace;
  //   // exit(1);
  // }
  
  /*
   * Make array of already treated face 
   */
  int *AlreadyTreat = (int * ) malloc( sizeof(int *) * nFac);
  for(int i=0; i < nFac; i++){
    AlreadyTreat[i] = -1;
  }

  /*
   * Begin with the first face of the array
   */
  for(int iPos=0; iPos < nFac; iPos++){

    /** Get the current Idx in dFaceVtx **/
    int iFacIdx = dFaceVtxIdx[iAbsFace];

    /* 
     * Verbose 
     */
    // printf("---------------------------- : %d - %d - %d\n", iPos, IdxFace[iPos ], data[IdxFace[iPos ]]);
    
    /* 
     * Face deja traité ?? 
     */
    if(AlreadyTreat[iPos] != 1){

      int curFac = IdxFace[iPos ];
      int nVtx1  = data[curFac+2];

      /** Prepare ElemCon **/
      int *ElmCon1 = (int *) malloc( sizeof(int *) * nVtx1); /* To sort out of loop */
      int  tKey1   = 0;
      for(int iVtx=0; iVtx < nVtx1; iVtx++){
        ElmCon1[iVtx] = data[curFac+3+iVtx];
        tKey1        += data[curFac+3+iVtx];
      }
      _quickSort_int(ElmCon1, 0, nVtx1-1);

      if(0 == 1){
        for(int iVtx=0; iVtx<nVtx1; iVtx++){
          printf("ElmCon1[%d] : %d \n", iVtx, ElmCon1[iVtx]);
        }
      }

      /** Loop on candidates faces **/
      for(int iNex=iPos+1; iNex < nFac; iNex++){

        /* 
         * Verbose 
         */
        // printf("++++++++++++++++++++++++++++++++++ : %d - %d - %d \n", iNex, IdxFace[iNex ], data[IdxFace[iNex ]]);

        /* 
         * On sait que la prochaine a traiter est forcement plus grande 
         */
        int nexFac = IdxFace[iNex ]; 
        int nVtx2  = data[nexFac+2];        
        
        // printf("+++++++++++++++++++ : %d - %d \n", nVtx1, nVtx2);
        /*
         * First sort - Number of Vertex is different -> Si Vtx different c'est pas la meme
         */
        if(nVtx1 == nVtx2){

          /*
           * Allocate memory and copy 
           */
          int *ElmCon2 = (int *) malloc( sizeof(int *) * nVtx1);
          int  tKey2   = 0;
          for(int iVtx=0; iVtx < nVtx1; iVtx++){
            ElmCon2[iVtx] = data[nexFac+3+iVtx];
            tKey2        += data[nexFac+3+iVtx];
          }

          /*
           * Sort ElemCon2 
           */
          _quickSort_int(ElmCon2, 0, nVtx2-1);
         
          if(0 == 1){
            for(int iVtx=0; iVtx < nVtx1; iVtx++){
              printf("ElmCon2[%d] : %d \n", iVtx, ElmCon2[iVtx]);
            }
          }

          /** Compare **/
          int isSameFace = 1;
          for(int iVtx=0; iVtx < nVtx1; iVtx++){
             if(ElmCon1[iVtx] != ElmCon2[iVtx]){
                isSameFace = -1;
                break;
             }
          }
          
          /* 
           * Verbose 
           */
          if(0 == 1){
            printf("tKey : %d / %d \n", tKey1, tKey2);
            assert(tKey1 == tKey2);
          }

          /*
           * Fill the Face data if Face is the same 
           */
          if(isSameFace == 1){

            /*
             * Stockage de parent elements 
             *    -> TODOUX -> Manage sign if Entrante /Sortante
             *    Rendre un contenu explicit : - si sortant // + si entrant 
             *    Il faut que le parent element d'une boundary soit le droit ...
             */
            // if(data[curFac] == -1){
            // if(data[curFac] > nCell ){
            //   dFaceCell[2*iAbsFace  ] = data[nexFac];
            //   // dFaceCell[2*iAbsFace+1] = 0;//data[curFac];
            //   dFaceCell[2*iAbsFace+1] = data[curFac];
  
            //   /* 
            //    * Fill FaceVtx connectivity array 
            //    */
            //   for(int iVtx=0; iVtx < nVtx1; iVtx++){
            //     dFaceVtx[iFacIdx+iVtx] = data[nexFac+3+iVtx];
            //   }
            //   printf("Impossible case now \n");
            //   exit(1);
            // }
            // else 
            // {
              dFaceCell[2*iAbsFace  ] = data[curFac];
              dFaceCell[2*iAbsFace+1] = data[nexFac];
  
              /* 
               * Fill FaceVtx connectivity array 
               */
              for(int iVtx=0; iVtx < nVtx1; iVtx++){
                dFaceVtx[iFacIdx+iVtx] = data[curFac+3+iVtx];
              }
            // }

            /*
             * Stockage FaceVtx connectivity Index
             */
            // printf("dFaceVtxIdx est faux %d \n", nVtx1);
            // dFaceVtxIdx[++iAbsFace] = nVtx1-1;
            // printf("dFaceVtxIdx est faux [%d] -> %d \n", iAbsFace+1, dFaceVtxIdx[iAbsFace+1]);
            dFaceVtxIdx[iAbsFace+1] = dFaceVtxIdx[iAbsFace] + nVtx1;
            iAbsFace++;

            /* 
             * Flags the two faces as treated 
             */
            AlreadyTreat[iPos] = 1;
            AlreadyTreat[iNex] = 1;

          } /** End if (isSameFace) **/

          /** Free Cendidate ElemCon **/
          free(ElmCon2);

        } /* Enf if Same Number of Vertex */

      } 
      /** Free current ElemCon **/
      free(ElmCon1);

    } /** End If alreadyTreat **/

    /* Boundary management **/
    if(AlreadyTreat[iPos] != 1){
      
      // printf("iAbsFace : %d \n", iAbsFace);
      // printf("iPos : %d \n", iPos);
      // printf("----------: %d - %d - %d\n", curFac, nVtx1, iFacIdx);

      int iFacIdx2 = dFaceVtxIdx[iAbsFace];
      int curFac  = IdxFace[iPos];
      int nVtx1   = data[curFac+2];
    
      dFaceCell[2*iAbsFace  ] = data[curFac];
      dFaceCell[2*iAbsFace+1] = 0;
    
      for(int iVtx=0; iVtx < nVtx1; iVtx++){
        dFaceVtx[iFacIdx2+iVtx] = data[curFac+3+iVtx];
      }
      
      /* Update index */
      dFaceVtxIdx[iAbsFace+1] = dFaceVtxIdx[iAbsFace] + nVtx1;
      iAbsFace++;
    
    }

  }

  /** Boundary management - Temporary **/


  /** Free **/
  free(AlreadyTreat);

  return iAbsFace;
}

static void
_make_absolute_face_numbering(PDM_DMesh_nodal_t* mesh)
{
  
  PDM_g_num_t nFaceProc = mesh->dNFace;
  PDM_g_num_t beg_NumAbs;
  
  PDM_MPI_Scan(&nFaceProc, &beg_NumAbs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, mesh->pdm_mpi_comm);
  beg_NumAbs -= nFaceProc;
  
  
  /** Compute the distribution of elements amont proc **/
  mesh->face_distrib = (PDM_g_num_t *) malloc((mesh->n_proc+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t _dNFace = (PDM_g_num_t) mesh->dNFace;
  PDM_MPI_Allgather((void *) &_dNFace,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&mesh->face_distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    mesh->pdm_mpi_comm);

  // mesh->face_distrib[0] = 1;
  mesh->face_distrib[0] = 0;

  for (int i = 1; i < mesh->n_proc+1; i++) {
    mesh->face_distrib[i] +=  mesh->face_distrib[i-1];
  }

  if (0 == 1) {
    printf("beg_NumAbs::Face : "PDM_FMT_G_NUM" \n", beg_NumAbs);
    printf("mesh->face_distrib : "PDM_FMT_G_NUM,  mesh->face_distrib[0]);
    for (int i = 1; i < mesh->n_proc+1; i++) {
      printf(" "PDM_FMT_G_NUM, mesh->face_distrib[i]);
    }
    printf("\n");
  }
}

/**
 * 
 * \brief Free vtx structure
 *
 * \param[inout]  vtx    Vertices
 *
 * \return        NULL
 * 
 */

static 
PDM_DMesh_nodal_vtx_t *
_vtx_free
(
 PDM_DMesh_nodal_vtx_t *vtx
)
{
  
  if (vtx != NULL) {

    if (vtx->distrib != NULL) {
      free (vtx->distrib);
      vtx->distrib = NULL;
    }

    free (vtx);
  }
  return NULL;
}


/**
 * 
 * \brief Initialize a mesh
 *
 * \param [inout]  mesh        Mesh
 * \param [in]     n_part      Number of partitions
 */

static void   
_mesh_init
(
PDM_DMesh_nodal_t *mesh,
const PDM_MPI_Comm comm,     
      PDM_g_num_t  nVtx,         
      PDM_g_num_t  nCel         
)
{
  int n_proc;
  int i_proc;
  
  PDM_MPI_Comm_size (comm, &n_proc);
  PDM_MPI_Comm_rank (comm, &i_proc);

  mesh->n_proc                   = n_proc;
  mesh->i_proc                   = i_proc;

  mesh->n_som_abs                = nVtx;          
  mesh->n_cell_abs               = nCel;          

  mesh->vtx                      = malloc(sizeof(PDM_DMesh_nodal_vtx_t ));
  mesh->vtx->_coords             = NULL;
  mesh->vtx->n_vtx               = 0;

  mesh->sections_std             = NULL;  
  mesh->sections_poly2d          = NULL;            
  mesh->sections_poly3d          = NULL;           

  mesh->pdm_mpi_comm             = comm;
  mesh->sections_id              = NULL;
  mesh->n_sections               = 0;

  mesh->n_dcell                  = -1;
  mesh->dCellFace                = NULL;
  mesh->dCellFaceIdx             = NULL;
  mesh->cell_distrib             = NULL;

  mesh->dNFace                   = -1;
  mesh->_dFaceVtx                = NULL;
  mesh->_dFaceVtxIdx             = NULL;
  mesh->_dFaceCell               = NULL;
  mesh->face_distrib             = NULL;
  
}

/**
 * 
 * \brief Update sections identifier list
 *
 * \param [inout]  mesh        Mesh
 */

static void   
_update_sections_id
(
PDM_DMesh_nodal_t *mesh
)
{
  int n_sections = 0;
  
  if (mesh->sections_std != NULL) {
    n_sections += PDM_Handles_n_get (mesh->sections_std);
  }
    
  if (mesh->sections_poly2d != NULL) {
    n_sections += PDM_Handles_n_get (mesh->sections_poly2d);
  }

  if (mesh->sections_poly3d != NULL) {
    n_sections += PDM_Handles_n_get (mesh->sections_poly3d);
  }
  
  if (mesh->n_sections < n_sections) {
    mesh->sections_id = (int *) realloc(mesh->sections_id, sizeof(int) * n_sections);
  }
  
  int k = 0;
  if (mesh->sections_std != NULL) {
    const int *id1 = PDM_Handles_idx_get (mesh->sections_std);
    int n = PDM_Handles_n_get (mesh->sections_std);
    for (int i = 0; i < n; i++) {
      mesh->sections_id[k++] = id1[i] + PDM_BLOCK_ID_BLOCK_STD;
    }
  }
    
  if (mesh->sections_poly2d != NULL) {
    const int *id1 = PDM_Handles_idx_get (mesh->sections_poly2d);
    int n = PDM_Handles_n_get (mesh->sections_poly2d);
    for (int i = 0; i < n; i++) {
      mesh->sections_id[k++] = id1[i] + PDM_BLOCK_ID_BLOCK_POLY2D;
    }
  }
  
  if (mesh->sections_poly3d != NULL) {
    const int *id1 = PDM_Handles_idx_get (mesh->sections_poly3d);
    int n = PDM_Handles_n_get (mesh->sections_poly3d);
    for (int i = 0; i < n; i++) {
      mesh->sections_id[k++] = id1[i] + PDM_BLOCK_ID_BLOCK_POLY3D;
    }
  }
  
  mesh->n_sections = n_sections;  

}


/**
 * 
 * \brief Free a standard section
 *
 * \param [inout]  _bloc_std    Standard section
 *
 * \return         Null
 *   
 */

static
PDM_DMesh_nodal_section_std_t *
_section_std_free
(
PDM_DMesh_nodal_section_std_t *_section_std
)
{

  if (_section_std == NULL) {
    return NULL;
  }
  
  if (_section_std->distrib != NULL) {
    free (_section_std->distrib);
    _section_std->distrib = NULL;
  }
 
  free(_section_std);
  return NULL;
}


/**
 * 
 * \brief Free a polygon section
 *
 * \param [inout]  _bloc_poly2d    Polygon section
 *
 * \return         Null
 *   
 */

static
PDM_DMesh_nodal_section_poly2d_t *
_section_poly2d_free
(
PDM_DMesh_nodal_section_poly2d_t *_section_poly2d
)
{
  
  if (_section_poly2d == NULL) {
    return NULL;
  }
  
  if (_section_poly2d->distrib != NULL) {
    free (_section_poly2d->distrib);
    _section_poly2d->distrib = NULL;
  }
 
  
  free(_section_poly2d);

  return NULL;
}


/**
 * 
 * \brief Free a polyhedron section
 *
 * \param [inout]  _section_poly3d    Polyhedron section
 *
 * \return         Null
 *   
 */

static
PDM_DMesh_nodal_section_poly3d_t *
_section_poly3d_free
(
PDM_DMesh_nodal_section_poly3d_t *_section_poly3d
)
{
  if (_section_poly3d == NULL) {
    return NULL;
  }
  
  if (_section_poly3d->distrib != NULL) {
    free (_section_poly3d->distrib);
    _section_poly3d->distrib = NULL;
  }
 
  
  free(_section_poly3d);

  return NULL;

}



///**
// *
// * \brief _section_elt_faces_get
// *
// * \param [in]     mesh               Current mesh
// * \param [in]     id_section         Section identifier
// * \param [inout]  elt_face_vtx_idx   Index of element faces connectivity (preallocated)
// * \param [inout]  elt_face_vtx       Element faces connectivity (preallocated)
// *
// */

//static int
//_section_size_elt_faces_get
//(
//      PDM_DMesh_nodal_t *mesh,
//      int               *s_elt_face_vtx_idx,
//      int               *s_elt_face_vtx,
//      int               *s_elt_face_cell
//)
//{
  
//  int _s_elt_face_vtx_idx = 0;
//  int _s_elt_face_vtx = 0;
  
  
//  int n_sections_std = PDM_Handles_n_get (mesh->sections_std);
//  const int *list_ind = PDM_Handles_idx_get (mesh->sections_std);
  
//  for (int i = 0; i < n_sections_std; i++) {
//    PDM_DMesh_nodal_section_std_t *section =
//      (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, list_ind[i]);
//    int n_face_elt = 0;
//    int n_sum_vtx_face = 0;
    
//    switch (section->t_elt) {
//    case PDM_MESH_NODAL_TRIA3:
//      n_face_elt = 3;
//      n_sum_vtx_face = 6;
//      break;
//    case PDM_MESH_NODAL_TETRA4:
//      n_face_elt = 4;
//      n_sum_vtx_face = 12;
//      break;
//    case PDM_MESH_NODAL_QUAD4:
//      n_face_elt = 4;
//      n_sum_vtx_face = 8;
//      break;
//    case PDM_MESH_NODAL_HEXA8:
//      n_face_elt = 6;
//      n_sum_vtx_face = 24;
//      break;
//    case PDM_MESH_NODAL_PYRAMID5:
//      n_face_elt = 5;
//      n_sum_vtx_face = 16;
//      break;
//    case PDM_MESH_NODAL_PRISM6:
//      n_face_elt = 5;
//      n_sum_vtx_face = 18;
//      break;
//    default:
//      PDM_error(__FILE__, __LINE__, 0, "Error _section_size_elt_faces_get : Element type is not taking int account\n");
//      exit(EXIT_FAILURE);
//    }
    
//    _s_elt_face_vtx_idx += section->n_elt * n_face_elt;
//    _s_elt_face_vtx += section->n_elt * n_sum_vtx_face;
//  }
      
//  int n_sections_poly2d = PDM_Handles_n_get (mesh->sections_poly2d);
//  list_ind = PDM_Handles_idx_get (mesh->sections_poly2d);
  
//  for (int i = 0; i < n_sections_poly2d; i++) {
//    PDM_DMesh_nodal_section_poly2d_t *section =
//      (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, list_ind[i]);
//    _s_elt_face_vtx_idx += section->_connec_idx[section->n_elt];
//    _s_elt_face_vtx += 2 * section->_connec_idx[section->n_elt];
//  }
      
//  int n_sections_poly3d = PDM_Handles_n_get (mesh->sections_poly3d);
//  list_ind = PDM_Handles_idx_get (mesh->sections_poly3d);
  
//  for (int i = 0; i < n_sections_poly3d; i++) {
//    PDM_DMesh_nodal_section_poly3d_t *section =
//      (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, list_ind[i]);
//    _s_elt_face_vtx_idx += section->n_face;
//    _s_elt_face_vtx += section->_facvtx[section->_facvtx_idx[section->n_face]];
//  }

//  *s_elt_face_cell =  _s_elt_face_vtx_idx;
//  *s_elt_face_vtx_idx = _s_elt_face_vtx_idx + 1;
//  *s_elt_face_vtx = _s_elt_face_vtx + 1;
  

//  return *s_elt_face_vtx - 1;
  
//}


///**
// *
// * \brief _section_elt_faces_get
// *
// * \param [in]     mesh               Current mesh
// * \param [in]     id_section         Section identifier
// * \param [inout]  elt_face_vtx_idx   Index of element faces connectivity (preallocated)
// * \param [inout]  elt_face_vtx       Element faces connectivity (preallocated)
// *
// */

//static void
//_section_elt_faces_add
//(
//      PDM_DMesh_nodal_t *mesh,
//const int                id_section,
//      int               *n_elt_current,
//      int               *n_face_current,
//      int               *elt_face_vtx_idx,
//      PDM_g_num_t       *elt_face_vtx,
//      PDM_g_num_t       *elt_face_cell
//)
//{
//  int _n_face_current = *n_face_current;
  
//  int         *_current_elt_face_vtx_idx = elt_face_vtx_idx + _n_face_current;
//  PDM_g_num_t *_current_elt_face_vtx     = elt_face_vtx + elt_face_vtx_idx[_n_face_current];
//  PDM_g_num_t *_current_elt_face_cell    = elt_face_cell +_n_face_current;
  
//  if (id_section < PDM_BLOCK_ID_BLOCK_POLY2D) {
  
//    PDM_DMesh_nodal_section_std_t *section =
//            (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, id_section);

//    if (section == NULL) {
//      PDM_error (__FILE__, __LINE__, 0, "Bad section identifier\n");
//    }

//    switch (section->t_elt) {

//    case PDM_MESH_NODAL_TRIA3:
//      {
//        const int n_face_elt        = 3;
//        const int n_sum_vtx_face    = 6;
//        const int n_sum_vtx_elt     = 3;

//        for (int ielt = 0; ielt < section->n_elt; ielt++) {
//          for (int iface = 0; iface < n_face_elt; iface++) {
//            _current_elt_face_vtx_idx[ielt * n_face_elt + iface + 1] =
//              _current_elt_face_vtx_idx[ielt * n_face_elt + iface] + 3;
//            _current_elt_face_cell[ielt * n_face_elt + iface    ] = *n_elt_current + ielt + 1;
//          }

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 0] = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 1] = section->_connec[n_sum_vtx_elt * ielt + 1];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 2] = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 3] = section->_connec[n_sum_vtx_elt * ielt + 2];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 4] = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 5] = section->_connec[n_sum_vtx_elt * ielt + 3];
//        }
        
//        *n_elt_current += section->n_elt;
//        *n_face_current += section->n_elt * n_face_elt;
//      }
//      break;

//    case PDM_MESH_NODAL_TETRA4:
//      {
//        const int n_face_elt        = 4;
//        const int n_sum_vtx_face    = 12;
//        const int n_sum_vtx_elt     = 4;

//        for (int ielt = 0; ielt < section->n_elt; ielt++) {
//          for (int iface = 0; iface < n_face_elt; iface++) {
//            _current_elt_face_vtx_idx[ielt * n_face_elt + iface + 1] =
//              _current_elt_face_vtx_idx[ielt * n_face_elt + iface] + 3;
//            _current_elt_face_cell[ielt * n_face_elt + iface    ] = *n_elt_current + ielt + 1;
//          }

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 0]  = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 1]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 2]  = section->_connec[n_sum_vtx_elt * ielt + 2];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 3]  = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 4]  = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 5]  = section->_connec[n_sum_vtx_elt * ielt + 1];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 6]  = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 7]  = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 8]  = section->_connec[n_sum_vtx_elt * ielt + 3];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 9]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 10] = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 11] = section->_connec[n_sum_vtx_elt * ielt + 2];

//        }
 
//        *n_elt_current += section->n_elt;
//        *n_face_current += section->n_elt * n_face_elt;

//    }
//      break;
//    case PDM_MESH_NODAL_QUAD4:
//      {
//        const int n_face_elt        = 4;
//        const int n_sum_vtx_face    = 8;
//        const int n_sum_vtx_elt     = 4;

//        for (int ielt = 0; ielt < section->n_elt; ielt++) {
//          for (int iface = 0; iface < n_face_elt; iface++) {
//            _current_elt_face_vtx_idx[ielt * n_face_elt + iface + 1] =
//              _current_elt_face_vtx_idx[ielt * n_face_elt + iface] + 3;
//            _current_elt_face_cell[ielt * n_face_elt + iface    ] = *n_elt_current + ielt + 1;
//          }

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 0] = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 1] = section->_connec[n_sum_vtx_elt * ielt + 1];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 2] = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 3] = section->_connec[n_sum_vtx_elt * ielt + 2];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 4] = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 5] = section->_connec[n_sum_vtx_elt * ielt + 3];
//        }
        
//        *n_elt_current  += section->n_elt;
//        *n_face_current += section->n_elt * n_face_elt;

//      }
//      break;
    
//    case PDM_MESH_NODAL_HEXA8:
//      {
//        const int n_face_elt        = 6;
//        const int n_sum_vtx_face    = 24;
//        const int n_sum_vtx_elt     = 8;

//        for (int ielt = 0; ielt < section->n_elt; ielt++) {
//          for (int iface = 0; iface < n_face_elt; iface++) {
//            _current_elt_face_vtx_idx[ielt * n_face_elt + iface + 1] =
//              _current_elt_face_vtx_idx[ielt * n_face_elt + iface] + 4;
//            _current_elt_face_cell[ielt * n_face_elt + iface    ] = *n_elt_current + ielt + 1;
//          }

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 0]  = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 1]  = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 2]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 3]  = section->_connec[n_sum_vtx_elt * ielt    ];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 4]  = section->_connec[n_sum_vtx_elt * ielt + 6];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 5]  = section->_connec[n_sum_vtx_elt * ielt + 7];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 6]  = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 7]  = section->_connec[n_sum_vtx_elt * ielt + 5];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 8]  = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 9]  = section->_connec[n_sum_vtx_elt * ielt + 7];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 10] = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 11] = section->_connec[n_sum_vtx_elt * ielt    ];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 12] = section->_connec[n_sum_vtx_elt * ielt + 7];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 13] = section->_connec[n_sum_vtx_elt * ielt + 6];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 14] = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 15] = section->_connec[n_sum_vtx_elt * ielt + 3];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 16] = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 17] = section->_connec[n_sum_vtx_elt * ielt + 6];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 18] = section->_connec[n_sum_vtx_elt * ielt + 5];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 19] = section->_connec[n_sum_vtx_elt * ielt + 1];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 20] = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 21] = section->_connec[n_sum_vtx_elt * ielt + 5];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 22] = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 23] = section->_connec[n_sum_vtx_elt * ielt + 0];
//        }
        
//        *n_elt_current  += section->n_elt;
//        *n_face_current += section->n_elt * n_face_elt;
        
//      }
//      break;
  
//    case PDM_MESH_NODAL_PYRAMID5:
//      {
//        const int n_face_elt        = 5;
//        const int n_sum_vtx_face    = 1*4 + 4*3;
//        const int n_sum_vtx_elt     = 5;

//        elt_face_vtx_idx[0] = 0;

//        for (int ielt = 0; ielt < section->n_elt; ielt++) {

//          for (int iface = 0; iface < n_face_elt; iface++) {
//            _current_elt_face_cell[ielt * n_face_elt + iface    ] = *n_elt_current + ielt + 1;
//          }

//          _current_elt_face_vtx_idx[ielt * n_face_elt + 1]  = elt_face_vtx_idx[ielt * n_face_elt    ] + 4;
//          _current_elt_face_vtx_idx[ielt * n_face_elt + 2]  = elt_face_vtx_idx[ielt * n_face_elt + 1] + 3;
//          _current_elt_face_vtx_idx[ielt * n_face_elt + 3]  = elt_face_vtx_idx[ielt * n_face_elt + 2] + 3;
//          _current_elt_face_vtx_idx[ielt * n_face_elt + 4]  = elt_face_vtx_idx[ielt * n_face_elt + 3] + 3;
//          _current_elt_face_vtx_idx[ielt * n_face_elt + 5]  = elt_face_vtx_idx[ielt * n_face_elt + 4] + 3;

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 0]  = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 1]  = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 2]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 3]  = section->_connec[n_sum_vtx_elt * ielt    ];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 4]  = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 5]  = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 6]  = section->_connec[n_sum_vtx_elt * ielt + 1];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 7]  = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 8]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 9]  = section->_connec[n_sum_vtx_elt * ielt + 2];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 10] = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 11] = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 12] = section->_connec[n_sum_vtx_elt * ielt + 3];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 13] = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 14] = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 15] = section->_connec[n_sum_vtx_elt * ielt    ];

//        }

//        *n_elt_current  += section->n_elt;
//        *n_face_current += section->n_elt * n_face_elt;

//      }
//      break;
    
//    case PDM_MESH_NODAL_PRISM6:
//      {
//        const int n_face_elt        = 5;
//        const int n_sum_vtx_face    = 3*4 + 2*3;
//        const int n_sum_vtx_elt     = 6;

//        for (int ielt = 0; ielt < section->n_elt; ielt++) {

//          for (int iface = 0; iface < n_face_elt; iface++) {
//            _current_elt_face_cell[ielt * n_face_elt + iface    ] = *n_elt_current + ielt + 1;
//          }

//          _current_elt_face_vtx_idx[ielt * n_face_elt + 1]  = elt_face_vtx_idx[ielt * n_face_elt    ] + 3;
//          _current_elt_face_vtx_idx[ielt * n_face_elt + 2]  = elt_face_vtx_idx[ielt * n_face_elt + 1] + 3;
//          _current_elt_face_vtx_idx[ielt * n_face_elt + 3]  = elt_face_vtx_idx[ielt * n_face_elt + 2] + 4;
//          _current_elt_face_vtx_idx[ielt * n_face_elt + 4]  = elt_face_vtx_idx[ielt * n_face_elt + 3] + 4;
//          _current_elt_face_vtx_idx[ielt * n_face_elt + 5]  = elt_face_vtx_idx[ielt * n_face_elt + 4] + 4;

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 0]  = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 1]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 2]  = section->_connec[n_sum_vtx_elt * ielt    ];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 3]  = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 4]  = section->_connec[n_sum_vtx_elt * ielt + 5];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 5]  = section->_connec[n_sum_vtx_elt * ielt + 3];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 6]  = section->_connec[n_sum_vtx_elt * ielt + 5];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 7]  = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 8]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 9]  = section->_connec[n_sum_vtx_elt * ielt + 2];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 10] = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 11] = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 12] = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 13] = section->_connec[n_sum_vtx_elt * ielt + 1];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 14] = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 15] = section->_connec[n_sum_vtx_elt * ielt + 5];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 16] = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 17] = section->_connec[n_sum_vtx_elt * ielt    ];

//        }

//        *n_elt_current  += section->n_elt;
//        *n_face_current += section->n_elt * n_face_elt;
    
//      }
//      break;
//    }
//  }
  
//  else if (id_section < PDM_BLOCK_ID_BLOCK_POLY3D) {

//    int _id = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
//    PDM_DMesh_nodal_section_poly2d_t *section =
//            (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, _id);
    
//    if (section == NULL) {
//      PDM_error (__FILE__, __LINE__, 0, "Bad section identifier\n");
//    }

//    int idx = 0;

//    for (int ielt = 0; ielt < section->n_elt; ielt++) {
//      int n_face_elt = section->_connec_idx[section->n_elt];
//      *n_face_current += n_face_elt;
//      int idx2 = section->_connec_idx[ielt];
//      for (int iface = 0; iface < n_face_elt; iface++) {
//        _current_elt_face_vtx_idx[idx + 1]  = _current_elt_face_vtx_idx[idx] + 2;
//        int inext = (iface + 1) % n_face_elt;
//        _current_elt_face_vtx[2 * idx    ]  = section->_connec[idx2 + iface];
//        _current_elt_face_vtx[2 * idx + 1]  = section->_connec[idx2 + inext];
//        _current_elt_face_cell[idx   ]  = *n_elt_current + ielt + 1;
//        idx += 1;
//      }
//    }
    
//    *n_elt_current += section->n_elt;
    
//  }
  
//  else {
    
//    PDM_error (__FILE__, __LINE__, 0, "PDM_BLOCK_ID_BLOCK_POLY3D : Not implemented yet\n");
 
//    //TODO: Compliqué car il faut faire des échanges Block_to_part pour recuperer les
//    // definitions des faces
//    // Il faut redupliquer les faces et les stocker comme pour les autres types de
//    // section
//    // Attention : Dans le cas d'un maillage avec une seule section poly3d, il ne faut
//    // rien faire.et ne pas passer dans cette fonction
    
//  }
  
//}
  

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   n_part   Number of partition on the current process
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */

int 
PDM_DMesh_nodal_create
(
const PDM_MPI_Comm comm,     
      PDM_g_num_t  nVtx,     
      PDM_g_num_t  nCel     
)
{
  PDM_DMesh_nodal_t *mesh = (PDM_DMesh_nodal_t *) malloc (sizeof(PDM_DMesh_nodal_t));

   if (sizeof(int) != sizeof(PDM_g_num_t)) {

    
    printf("PDM_DMesh_nodal : Erreur : Cette fontion ne fonctionne pas en 64bit\n");
    exit(1);
    
  }
  
  _mesh_init (mesh, comm, nVtx, nCel);
  
  if (mesh_handles == NULL) {
    mesh_handles = PDM_Handles_create (4);
  }
  
  return PDM_Handles_store (mesh_handles, mesh);
}


/**
 * \brief Free a nodal mesh structure
 *
 * \param [in]  idx   Nodal mesh handle
 *
 * \return      NULL
 *
 */

void
PDM_DMesh_nodal_free
(
const int hdl
)
{
  
  PDM_DMesh_nodal_t * mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh != NULL) {

    if (mesh->sections_id != NULL) {
      free (mesh->sections_id);
    }
    
    mesh->sections_id = NULL;
    
     _vtx_free (mesh->vtx);

    /* free standard sections */

    if (mesh->sections_std != NULL) {
      int n_sections_std = PDM_Handles_n_get (mesh->sections_std);
      const int *list_ind = PDM_Handles_idx_get (mesh->sections_std);

      while (n_sections_std > 0) {
        PDM_DMesh_nodal_section_std_t *_bloc_std = 
          (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, list_ind[0]);
        _section_std_free(_bloc_std);
        PDM_Handles_handle_free (mesh->sections_std, list_ind[0], PDM_FALSE);
        n_sections_std = PDM_Handles_n_get (mesh->sections_std);
      }

      mesh->sections_std = PDM_Handles_free (mesh->sections_std); 
    }

    /* Free polygon sections */ 

    if (mesh->sections_poly2d != NULL) {
      int n_sections_poly2d = PDM_Handles_n_get (mesh->sections_poly2d);
      const int *list_ind = PDM_Handles_idx_get (mesh->sections_poly2d);

      while (n_sections_poly2d > 0) {
        PDM_DMesh_nodal_section_poly2d_t *_bloc_poly2d = 
          (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, list_ind[0]);
        _section_poly2d_free(_bloc_poly2d);
        PDM_Handles_handle_free (mesh->sections_poly2d, list_ind[0], PDM_FALSE);
        n_sections_poly2d = PDM_Handles_n_get (mesh->sections_poly2d);
      }

      mesh->sections_poly2d = PDM_Handles_free (mesh->sections_poly2d); 
    }

    /* Free polyhedron sections */ 

    if (mesh->sections_poly3d != NULL) {
      int n_sections_poly3d = PDM_Handles_n_get (mesh->sections_poly3d);
      const int *list_ind = PDM_Handles_idx_get (mesh->sections_poly3d);

      while (n_sections_poly3d > 0) {
        PDM_DMesh_nodal_section_poly3d_t *_bloc_poly3d = 
          (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, list_ind[0]);
        _section_poly3d_free(_bloc_poly3d);
        PDM_Handles_handle_free (mesh->sections_poly3d, list_ind[0], PDM_FALSE);
        n_sections_poly3d = PDM_Handles_n_get (mesh->sections_poly3d);
      }

      mesh->sections_poly3d = PDM_Handles_free (mesh->sections_poly3d); 
    }
    
    if (mesh->sections_id != NULL) {
      free (mesh->sections_id);
    }

    if (mesh->dCellFaceIdx != NULL) {
      free (mesh->dCellFaceIdx);
    }
     
    if (mesh->dCellFace != NULL) {
      free (mesh->dCellFace);
    }
     
    if (mesh->cell_distrib != NULL) {
      free (mesh->cell_distrib);
    }
     
    if (mesh->_dFaceVtxIdx != NULL) {
      free (mesh->_dFaceVtxIdx);
    }
     
    if (mesh->_dFaceVtx != NULL) {
      free (mesh->_dFaceVtx);
    }
     
    if (mesh->face_distrib != NULL) {
      free (mesh->face_distrib);
    }

     free(mesh);

    PDM_Handles_handle_free (mesh_handles, hdl, PDM_FALSE);
  
    int n_mesh_array = PDM_Handles_n_get (mesh_handles); 

    if (n_mesh_array == 0) {
      mesh_handles = PDM_Handles_free (mesh_handles);
    }
  }
}


/**
 * \brief Define partition vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 * \param [in]  n_vtx     Number of vertices
 * \param [in]  coords    Interlaced coordinates (size = 3 * \ref n_vtx)
 *
 */

void
PDM_DMesh_nodal_coord_set
(
 const int          hdl,
 const int          n_vtx,  
 const PDM_real_t  *coords  
)
{
  PDM_DMesh_nodal_t * mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  PDM_DMesh_nodal_vtx_t *vtx = mesh->vtx;

  if (vtx->_coords != NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Vertices are already defined\n");
  }

  /* Mapping memoire */

  vtx->n_vtx   = n_vtx;
  vtx->_coords = coords;

  vtx->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mesh->n_proc + 1));
  
  PDM_g_num_t *_distrib = vtx->distrib + 1;
  _distrib[0] = 0;
  PDM_g_num_t _n_vtx = n_vtx;
  
  PDM_MPI_Scan (&_n_vtx, _distrib, 1, PDM__PDM_MPI_G_NUM, 
                PDM_MPI_SUM, mesh->pdm_mpi_comm);


}


/**
 * \brief  Return number of vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Number of vertices
 *
 */

int
PDM_DMesh_nodal_n_vtx_get
(
 const int          hdl
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  PDM_DMesh_nodal_vtx_t *vtx = mesh->vtx;

  return vtx->n_vtx;
}


/**
 * \brief  Return coordinates of vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Coordinates of vertices
 *
 */

const double *
PDM_DMesh_nodal_vtx_get
(
 const int          hdl
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  PDM_DMesh_nodal_vtx_t *vtx = mesh->vtx;

  return vtx->_coords;
}


/**
 * \brief  Return number of sections
 *
 * \param [in]  hdl            Nodal mesh handle
 *
 * \return  Number of sections
 *
 */

int
PDM_DMesh_nodal_n_sections_get
(
 const int   hdl
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  
  return mesh->n_sections;

}


/**
 * \brief  Return sections identifier
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 *
 * \return  Blocks identifier
 *
 */

int *
PDM_DMesh_nodal_sections_id_get
(
const int   hdl
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  return mesh->sections_id;

}


/**
 * \brief  Return type of section
 *
 * \param [in]  hdl        Distributed nodal mesh handle
 * \param [in]  id_section   Block identifier
 *
 * \return  Type of section
 *
 */

PDM_Mesh_nodal_elt_t
PDM_DMesh_nodal_section_type_get
(
const int   hdl,
const int   id_section     
)
{
  PDM_DMesh_nodal_t *mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  PDM_Mesh_nodal_elt_t t_elt;
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  if (id_section < PDM_BLOCK_ID_BLOCK_POLY2D) {
  
    t_elt = PDM_MESH_NODAL_POLY_3D;
    PDM_DMesh_nodal_section_std_t *section = 
            (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, id_section);

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad section identifier\n");  
    }

    t_elt = section->t_elt;
  }
  
  else if (id_section < PDM_BLOCK_ID_BLOCK_POLY3D) {
    
    t_elt = PDM_MESH_NODAL_POLY_2D;

  }
  
  else {
    
    t_elt = PDM_MESH_NODAL_POLY_3D;

  }
  
  return t_elt;
    
}


/**
 * \brief  Add a new section to the current mesh
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  st_free_data   Status of Release of the memory 
 *                             when the section is destroyed
 * \param [in]  id_section       Block identifier
 *
 * \return Block identifier     
 *
 */

int 
PDM_DMesh_nodal_section_add 
(
const int                    hdl,
const PDM_Mesh_nodal_elt_t   t_elt
)
{

  PDM_DMesh_nodal_t *mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  int id_section;
  
  switch (t_elt) {

  case PDM_MESH_NODAL_POINT    :    
  case PDM_MESH_NODAL_BAR2     :   
  case PDM_MESH_NODAL_TRIA3    :    
  case PDM_MESH_NODAL_QUAD4    :    
  case PDM_MESH_NODAL_TETRA4   :     
  case PDM_MESH_NODAL_PYRAMID5 :     
  case PDM_MESH_NODAL_PRISM6   :     
  case PDM_MESH_NODAL_HEXA8    : 
    {
      /* Mise a jour du tableau de stockage */

      if (mesh->sections_std == NULL) {
        mesh->sections_std = PDM_Handles_create (4);
      } 
      
      /* Allocation du bloc */
      
      PDM_DMesh_nodal_section_std_t *section_std = (PDM_DMesh_nodal_section_std_t *) malloc(sizeof(PDM_DMesh_nodal_section_std_t));

      id_section = PDM_Handles_store (mesh->sections_std, section_std);

      /* Intialisation du bloc */

      section_std->t_elt = t_elt;

      section_std->n_elt      = -1;
      section_std->_connec    = NULL;
      section_std->distrib    = NULL;
      
      id_section += PDM_BLOCK_ID_BLOCK_STD;
      if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {
        PDM_error(__FILE__, __LINE__, 0, "The number of standard sections must be less than %d\n", 
               PDM_BLOCK_ID_BLOCK_POLY2D);
        abort();
      }
    }
    
    break;

  case PDM_MESH_NODAL_POLY_2D  :    
    {
      /* Mise a jour du tableau de stockage */

      if (mesh->sections_poly2d == NULL) {
        mesh->sections_poly2d = PDM_Handles_create (4);
      }
      
      /* Allocation du bloc */
      
      PDM_DMesh_nodal_section_poly2d_t *section_poly2d =
              (PDM_DMesh_nodal_section_poly2d_t *) malloc(sizeof(PDM_DMesh_nodal_section_poly2d_t));

      id_section = PDM_Handles_store (mesh->sections_poly2d, section_poly2d);
      
      /* Intialisation du bloc */

      section_poly2d->n_elt       = -1;
      section_poly2d->_connec_idx = NULL;
      section_poly2d->_connec     = NULL;
      section_poly2d->distrib     = NULL;

      id_section += PDM_BLOCK_ID_BLOCK_POLY2D;
      if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {
        PDM_error(__FILE__, __LINE__, 0, "The number of polygon sections must be less than %d\n",
               PDM_BLOCK_ID_BLOCK_POLY3D - PDM_BLOCK_ID_BLOCK_POLY2D);
      }
    }

    break;

  case PDM_MESH_NODAL_POLY_3D  :     
    {
      /* Mise a jour du tableau de stockage */

      if (mesh->sections_poly3d == NULL) {
        mesh->sections_poly3d = PDM_Handles_create (4);
      }
      
      /* Allocation du bloc */
      
      PDM_DMesh_nodal_section_poly3d_t *section_poly3d =
              (PDM_DMesh_nodal_section_poly3d_t *) malloc(sizeof(PDM_DMesh_nodal_section_poly3d_t));

      id_section = PDM_Handles_store (mesh->sections_poly3d, section_poly3d);

      /* Intialisation du bloc */

      section_poly3d->n_elt        = -1;
      section_poly3d->n_face       = -1;
      section_poly3d->_facvtx_idx  = NULL;
      section_poly3d->_facvtx      = NULL;
      section_poly3d->_cellfac_idx = NULL;
      section_poly3d->_cellfac     = NULL;
      section_poly3d->distrib      = NULL;

      id_section += PDM_BLOCK_ID_BLOCK_POLY3D;

    }

    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown element type\n");

  }

  _update_sections_id (mesh);
  return id_section ;
  
}


/**
 * \brief Define a standard section
 *
 *  - PDM_MESH_NODAL_POINT :
 *
 *   1 x            
 *
 *  - PDM_MESH_NODAL_BAR2 :
 *
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_TRIA3 :   
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - PDM_MESH_NODAL_QUAD4 :          
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - PDM_MESH_NODAL_TETRA4 :    
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - PDM_MESH_NODAL_PYRAMID5 :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_PRSIM6 :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  - PDM_MESH_NODAL_HEXA8 :   
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect        Connectivity
 *
 */

void
PDM_DMesh_nodal_section_std_set 
(
const int          hdl,
const int          id_section,     
const int          n_elt,    
      PDM_g_num_t *connec   
)
{
  PDM_DMesh_nodal_t *mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
  
  PDM_printf(" PDM_Mesh_nodal_section_std_set - _id_section : %i  \n ", _id_section);
  PDM_printf(" PDM_Mesh_nodal_section_std_set - n_elt       : %i  \n ", n_elt);
  
  PDM_DMesh_nodal_section_std_t *section = (PDM_DMesh_nodal_section_std_t *) 
     PDM_Handles_get (mesh->sections_std, _id_section);
  
  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }
   
  /* Mapping */

  section->n_elt   = n_elt;
  section->_connec = connec;

  section->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mesh->n_proc + 1));
  
  /* Creation of distribution */
  PDM_g_num_t begNumAbs;
  PDM_g_num_t _n_elt = n_elt;
  
  PDM_MPI_Scan (&_n_elt, &begNumAbs, 1, PDM__PDM_MPI_G_NUM, 
                PDM_MPI_SUM, mesh->pdm_mpi_comm);
  begNumAbs -= _n_elt;

  PDM_MPI_Allgather((void *) &_n_elt,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&section->distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    mesh->pdm_mpi_comm);
  
  section->distrib[0] = 0;
  for (int i = 1; i < mesh->n_proc + 1; i++) {
    section->distrib[i] +=  section->distrib[i-1];
  }
  
}


/**
 * \brief Return standard section description
 *
 *  - PDM_MESH_NODAL_POINT :
 *
 *   1 x            
 *
 *  - PDM_MESH_NODAL_BAR2 :
 *
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_TRIA3 :   
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - PDM_MESH_NODAL_QUAD4 :          
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - PDM_MESH_NODAL_TETRA4 :    
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - PDM_MESH_NODAL_PYRAMID5 :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_PRSIM6 :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  - PDM_MESH_NODAL_HEXA8 :   
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 *
 * \return  connect           Connectivity
 *
 */

PDM_g_num_t *   
PDM_DMesh_nodal_section_std_get 
(   
const int            hdl,
const int            id_section     
) 
{
  PDM_DMesh_nodal_t *mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
  
  const PDM_DMesh_nodal_section_std_t *section = 
  (const PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, _id_section);
  
  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }
  
  return section->_connec;
}


/**
 * \brief Get number of section elements
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  id_part        Partition identifier
 *
 * \return      Number of elements
 *  
 */


/**
 * \brief Get number of section elements
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 *
 * \return      Number of elements
 *  
 */

int
PDM_DMesh_nodal_section_n_elt_get 
(   
const int            hdl,
const int            id_section     
) 
{
  PDM_DMesh_nodal_t *mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  int _id_section;
  
  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {
  
    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;
  
    const PDM_DMesh_nodal_section_poly3d_t *section = 
    (const PDM_DMesh_nodal_section_poly3d_t *) 
     PDM_Handles_get (mesh->sections_poly3d, _id_section);
  
    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
    }
      
    return section->n_elt;
  }
  
  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {
  
    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
  
    const PDM_DMesh_nodal_section_poly2d_t *section = 
    (const PDM_DMesh_nodal_section_poly2d_t *) 
     PDM_Handles_get (mesh->sections_poly2d, _id_section);
  
    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad polygon section identifier\n");
    }
  
    return section->n_elt;
  }
  
  else {
  
    _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
  
    const PDM_DMesh_nodal_section_std_t *section = 
    (const PDM_DMesh_nodal_section_std_t *) 
     PDM_Handles_get (mesh->sections_std, _id_section);
  
    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad polyhedron section identifier\n");
    }
  
    return section->n_elt;
  }
  
}


/**
 * \brief Define a polygon section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [in]  connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 *
 */
 
void
PDM_DMesh_nodal_section_poly2d_set 
(
const int            hdl,
const int            id_section, 
const PDM_l_num_t    n_elt,    
      PDM_l_num_t   *connec_idx,   
      PDM_g_num_t   *connec
)
{
  PDM_DMesh_nodal_t *mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
  
  PDM_DMesh_nodal_section_poly2d_t *section = 
          (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, _id_section);

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  /* Mapping */

  section->n_elt       = n_elt;
  section->_connec_idx = connec_idx;
  section->_connec     = connec;

  section->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mesh->n_proc + 1));
  
  /* Creation of distribution */
  PDM_g_num_t begNumAbs;
  PDM_g_num_t _n_elt = n_elt;
  
  PDM_MPI_Scan (&_n_elt, &begNumAbs, 1, PDM__PDM_MPI_G_NUM, 
                PDM_MPI_SUM, mesh->pdm_mpi_comm);
  begNumAbs -= _n_elt;

  PDM_MPI_Allgather((void *) &_n_elt,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&section->distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    mesh->pdm_mpi_comm);
  
  section->distrib[0] = 0;
  for (int i = 1; i < mesh->n_proc + 1; i++) {
    section->distrib[i] +=  section->distrib[i-1];
  }
  
}


/**
 * \brief Return a polygon section description
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [out] connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [out] connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 *
 */
 
void
PDM_DMesh_nodal_section_poly2d_get 
(
 const int          hdl,
 const int          id_section, 
       PDM_l_num_t  **connec_idx,   
       PDM_g_num_t  **connec
) 
{
  PDM_DMesh_nodal_t *mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
  
  PDM_DMesh_nodal_section_poly2d_t *section = 
          (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, _id_section);

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }
  
  *connec_idx = section->_connec_idx;
  *connec     = section->_connec;
 
}


/**
 * \brief Define a polyhedra section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_face         Number of faces used to describe polyhedra
 * \param [in]  facvtx_idx     Index of face vertex connectivity
 * \param [in]  facvtx         Face vertex connectivity
 * \param [in]  cellfac_idx    Index of cell face connectivity
 * \param [in]  cellfac        Cell face connectivity
 *
 */

void
PDM_DMesh_nodal_section_poly3d_set 
(
const int            hdl,
const int            id_section, 
const PDM_l_num_t    n_elt,    
const PDM_l_num_t    n_face,   
      PDM_l_num_t   *facvtx_idx,   
      PDM_g_num_t   *facvtx,
      PDM_l_num_t   *cellfac_idx,   
      PDM_g_num_t   *cellfac
)
{
  PDM_DMesh_nodal_t *mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
  
  PDM_DMesh_nodal_section_poly3d_t *section = 
          (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, _id_section);

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  section->n_elt        = n_elt;
  section->n_face       = n_face;
  section->_facvtx_idx  = facvtx_idx;
  section->_facvtx      = facvtx;
  section->_cellfac_idx = cellfac_idx;
  section->_cellfac     = cellfac;
  
  section->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mesh->n_proc + 1));
  
  /* Creation of distribution */
  PDM_g_num_t begNumAbs;
  PDM_g_num_t _n_elt = n_elt;
  
  PDM_MPI_Scan (&_n_elt, &begNumAbs, 1, PDM__PDM_MPI_G_NUM, 
                PDM_MPI_SUM, mesh->pdm_mpi_comm);
  begNumAbs -= _n_elt;

  PDM_MPI_Allgather((void *) &_n_elt,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&section->distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    mesh->pdm_mpi_comm);
  
  section->distrib[0] = 0;
  for (int i = 1; i < mesh->n_proc + 1; i++) {
    section->distrib[i] +=  section->distrib[i-1];
  }
  
  
}


/**
 * \brief Define a polyhedra section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [out]  n_face         Number of faces used to describe polyhedra
 * \param [out]  facvtx_idx     Index of face vertex connectivity
 * \param [out]  facvtx         Face vertex connectivity
 * \param [out]  cellfac_idx    Index of cell face connectivity
 * \param [out]  cellfac        Cell face connectivity
 *
 */

void
PDM_DMesh_nodal_section_poly3d_get 
(
const int            hdl,
const int            id_section, 
      PDM_l_num_t   *n_face,   
      PDM_l_num_t  **facvtx_idx,   
      PDM_g_num_t  **facvtx,
      PDM_l_num_t  **cellfac_idx,   
      PDM_g_num_t  **cellfac
)
{
  PDM_DMesh_nodal_t *mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
  
  PDM_DMesh_nodal_section_poly3d_t *section = 
          (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, _id_section);

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  *n_face      = section->n_face;
  *facvtx_idx  = section->_facvtx_idx;
  *facvtx      = section->_facvtx;
  *cellfac_idx = section->_cellfac_idx;
  *cellfac     = section->_cellfac;

}


/**
 * \brief  Return total number of elements of a distributed mesh
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Return number elements of a partition
 * 
 */

PDM_g_num_t
PDM_DMesh_nodal_total_n_cell_get
(
const int  hdl
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  PDM_g_num_t total_n_cell = 0;
  
  int n_sections_std = PDM_Handles_n_get (mesh->sections_std);
  const int *list_ind = PDM_Handles_idx_get (mesh->sections_std);
  
  for (int i = 0; i < n_sections_std; i++) {
    PDM_DMesh_nodal_section_std_t *_section_std = 
      (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, list_ind[i]);
    total_n_cell += _section_std->distrib[mesh->n_proc];
  }
      
  int n_sections_poly2d = PDM_Handles_n_get (mesh->sections_poly2d);
  list_ind = PDM_Handles_idx_get (mesh->sections_poly2d);
  
  for (int i = 0; i < n_sections_poly2d; i++) {
    PDM_DMesh_nodal_section_poly2d_t *_section_poly2d = 
      (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, list_ind[i]);
    total_n_cell += _section_poly2d->distrib[mesh->n_proc];
  }
      
  int n_sections_poly3d = PDM_Handles_n_get (mesh->sections_poly3d);
  list_ind = PDM_Handles_idx_get (mesh->sections_poly3d);
  
  for (int i = 0; i < n_sections_poly3d; i++) {
    PDM_DMesh_nodal_section_poly3d_t *_section_poly3d = 
      (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, list_ind[i]);
    total_n_cell += _section_poly3d->distrib[mesh->n_proc];
  }
  
  return total_n_cell;

}


/**
 * \brief  Return total number of faces of a distributed mesh
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Return total number of faces
 * 
 */

PDM_g_num_t
PDM_DMesh_nodal_total_n_face_get
(
const int  hdl
) 
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  // if (mesh->dCellFace == NULL) {
  //   PDM_error (__FILE__, __LINE__, 0, "Not implemented yet\n");
  // }
  
  return mesh->face_distrib[mesh->n_proc];
  
}


/**
 * \brief  Return total number of vertices of a distributed mesh
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Return total number of vertices
 * 
 */

PDM_g_num_t
PDM_DMesh_nodal_total_n_vtx_get
(
const int  hdl
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  PDM_DMesh_nodal_vtx_t *vtx = mesh->vtx;

  return vtx->distrib[mesh->n_proc];
}


/**
 * \brief  Compute cell->face connectivity
 *
 * \param [in]   hdl              Distributed nodal mesh handle
 * 
 */

void
PDM_DMesh_nodal_cell_face_compute
(
const int   hdl
)
{
  

  PDM_printf("PDM_DMesh_nodal_cell_face_compute \n ");
  
  /* Get current structure to treat */
  PDM_DMesh_nodal_t *mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  
  /* Verbose */
  if(1 == 0)
  {
    PDM_printf(" -------- \n ");
    PDM_printf(" n_proc     : %i \n ", mesh->n_proc);
    PDM_printf(" n_sections : %i \n ", mesh->n_sections);
    
    PDM_printf(" -------- \n ");
    
    // TODO : Sure ?? Check with Eric 
    
    for(int iSection=0; iSection < mesh->n_sections; iSection++){
      
      PDM_DMesh_nodal_section_std_t *sectionStd = (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, iSection);
      
      PDM_printf(" ------------------------------ \n ");
      PDM_printf(" iSection : %i  \n ", iSection);
      PDM_printf(" t_elt    : %i  \n ", sectionStd->t_elt);
      PDM_printf(" n_elt    : %i  \n ", sectionStd->n_elt);
      
      PDM_printf(" DitribElmt ... \n ");
      for(int iProc=0; iProc < mesh->n_proc + 1; iProc++)
        PDM_printf(PDM_FMT_G_NUM, sectionStd->distrib[iProc]);
      PDM_printf("\n");
    }
    
  }
  
  /* Creation of element distribution among all sections */
  mesh->sectionDistribution    = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mesh->n_sections + 1));
  mesh->sectionDistribution[0] = 0;
  
  for (int iSection = 0; iSection < mesh->n_sections; iSection++) {
    PDM_DMesh_nodal_section_std_t *sectionStd = (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, iSection);
    mesh->sectionDistribution[iSection+1] = sectionStd->distrib[mesh->n_proc];
  }
  for (int iSection = 1; iSection < mesh->n_sections + 1; iSection++) {
    mesh->sectionDistribution[iSection] +=  mesh->sectionDistribution[iSection-1];
  }
  
  /* Verbose */
  if(1 == 0)
  {
    PDM_printf(" ------------------------------ \n ");
    for(int iSection=0; iSection < mesh->n_sections+1; iSection++){
      PDM_printf("%i ", mesh->sectionDistribution[iSection]);
    }
  }
  
  
  /* Build an approximate number of faces for current processor */
  PDM_g_num_t dnElemTot  = 0;
  for (int iSection = 0; iSection < mesh->n_sections; iSection++) {    
    PDM_DMesh_nodal_section_std_t *sectionStd = (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, iSection);
    dnElemTot += ( sectionStd->distrib[mesh->i_proc+1] - sectionStd->distrib[mesh->i_proc] );
  }
  PDM_g_num_t nFacApprox = dnElemTot*8;
  
  if(1 == 0)
  {
    PDM_printf("nFacApprox : %i  \n ", nFacApprox);
  }
  
  
  /*
   * Allocate Disctribute Hash Key 
   */
  PDM_g_num_t* LNToGN = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * nFacApprox );
  
  
  /** Allocate (SurDim the part_data and part_stride ) **/
  int nDataApprox = dnElemTot*6*4*2;
  int* part_data = (int *) malloc( sizeof(int) * nDataApprox );
  int* part_stri = (int *) malloc( sizeof(int) * nFacApprox  );

  /** Initialisation of some data **/
  int nFace  = 0;
  int nData  = 0;
  
  /** Loop over Elements **/
  for (int iSection = 0; iSection < mesh->n_sections; iSection++) {        
    
    /* Get current section */
    PDM_DMesh_nodal_section_std_t *sectionStd = (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, iSection);
    
    /** Convert type to stride on elements **/
    int nVtxPerElmt = _get_size_of_element(sectionStd->t_elt);
    int nFacPerElmt = _get_nbface_per_element(sectionStd->t_elt);

    /** Get the index wich compute the elements **/
    int*  nVtxpFac  = (int * ) malloc( nFacPerElmt * sizeof(int * ) );
    int** tabFacVtx = (int **) malloc( nFacPerElmt * sizeof(int **) );

    /** Fill up nVtxpFac and tabFacVtx **/
    _get_elmt_info(sectionStd->t_elt, nVtxpFac, tabFacVtx);
    
    // int BegE = sectionStd->distrib[mesh->i_proc  ]*nVtxPerElmt  ; //+ offset;
    // int EndE = sectionStd->distrib[mesh->i_proc+1]*nVtxPerElmt-1; //+ offset;
      
    /** Verbose **/
    // if(0 == 0){
    //   printf("BegE/EndE   : %d/%d     \n", BegE, EndE);    
    // }
    
    /** Each section is composed of multiple elements **/
    for(int iElmt = sectionStd->distrib[mesh->i_proc]; iElmt < sectionStd->distrib[mesh->i_proc+1]; iElmt++){

      // int iOffSet = iElmt*nVtxPerElmt-sectionStd->distrib[mesh->i_proc];
      int iOffSet = (iElmt-sectionStd->distrib[mesh->i_proc])*nVtxPerElmt;
      // printf("[%i] - iOffset : %i --- %i \n",mesh->i_proc,  iOffSet, sectionStd->distrib[mesh->i_proc]);
      /** For each Elemt we decompose it to face  **/
      for(int iFace = 0; iFace < nFacPerElmt; iFace++){
        
        /* Compute the key */
        int iKey = _compute_key(sectionStd->_connec, tabFacVtx[iFace], iOffSet, nVtxpFac[iFace]);
        
        /* Build the LNToGN */
        // Attnetion ici c'est faux je pense !!!!
        if(nFace >= nFacApprox){
          nFacApprox = 2*nFacApprox;
          LNToGN    = (PDM_g_num_t *) realloc(LNToGN   , sizeof(PDM_g_num_t) * nFacApprox );
          part_stri = (int         *) realloc(part_stri, sizeof(int) * nFacApprox );
        }
        if(nData >= nDataApprox){
          nDataApprox = 2*nDataApprox;
          part_data = (int         *) realloc(part_data, sizeof(int) * nDataApprox );
        }
        
        LNToGN[nFace] = iKey;
        
        /** Prepare part_data - iElem iType nVtx (i1, i2, i3, ... ) - **/
        part_data[nData++] = iElmt+1+mesh->sectionDistribution[iSection];
        part_data[nData++] = 1;
        
        /* Setup nVtx */
        part_data[nData++] = nVtxpFac[iFace] ;
        
        for(int iVtx=0; iVtx < nVtxpFac[iFace]; iVtx++){
          // printf("lVtxIdx[%d] : %d -> %d \n ", iVtx, lVtxIdx[iVtx],dCellVtxArr[lVtxIdx[iVtx]]);
          // printf("lVtxIdx[%d] : %d -> %d \n ", BegE+tabFacVtx[iFace][iVtx], fetch->_dElemArr[BegE+tabFacVtx[iFace][iVtx]], iKey);
          // printf("lVtxIdx[%d] : %d -> %d \n ", iVtx, sectionStd->_connec[iOffSet+tabFacVtx[iFace][iVtx]], iKey);
          part_data[nData++] = sectionStd->_connec[iOffSet+tabFacVtx[iFace][iVtx]];
        }
  
        /** Prepare part_data **/
        part_stri[nFace] = 1+1+1+nVtxpFac[iFace];
  
        /** Go to next face **/
        nFace++;
      
      } /* End face loop */
    } /* End element loop */
    
    /** Free **/
    free(nVtxpFac);
    for(int i=0; i < nFacPerElmt; i++){
      free(tabFacVtx[i]);  
    }
    free(tabFacVtx);
    
  } /* End section loop */
  
  
  /*
   * Re-Allocate Array
   */
  LNToGN    = (PDM_g_num_t *) realloc((LNToGN   ), nFace * sizeof(PDM_g_num_t * ));
  part_stri = (int *) realloc((part_stri), nFace * sizeof(int * ));
  part_data = (int *) realloc((part_data), nData * sizeof(int * ));

  /* 
   * Verbose 
   */
  if(0 == 1){
    printf("------------------------------------------- \n");
    printf("nFace       : %d \n", nFace );
    printf("nData       : %d \n", nData );
    printf("nFacApprox  : %d \n", nFacApprox );
    printf("------------------------------------------- \n");
    for (int i = 0; i < nFace; i++) {
      printf("[%d] - i/iKey  : %d - %d \n", mesh->i_proc, i, LNToGN[i]);
    }
    printf("------------------------------------------- \n");
  }

  /* 
   * Create PartToBlock Structure 
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      &LNToGN,
                                                      &nFace,
                                                      1,
                                                      mesh->pdm_mpi_comm); 

  /* 
   * Exchange data in one call ( better do multiple -> Eric ??)
   */
  int *BlkStri = NULL;
  int *BlkData = NULL;
    
  int dataSize = PDM_part_to_block_exch(          ptb,
                                                  sizeof(PDM_g_num_t),
                                                  PDM_STRIDE_VAR,
                                                  -1,
                                                  &part_stri,
                                        (void **) &part_data,                                               
                                                  &BlkStri,
                                        (void **) &BlkData);

  /*
   *  Get the size of the current process bloc
   */
  int BlkSize = PDM_part_to_block_n_elt_block_get(ptb);

  /* 
   * Verbose 
   */
  if(0 == 1){
    printf("BlkSize : %d\n", BlkSize);
    for(int i = 0; i < BlkSize; i++) {
      // printf("BlockData[%d]    : %d \n", i, BlkData[i]);
      printf("BlkStri[%d] : %d \n", i, BlkStri[i]);
    }
  }
  
  /*
   *  Creation of array of diplacement
   */
  int* BlkStriIdx = (int *) malloc( sizeof(int *) * (BlkSize+1) ) ;
  BlkStriIdx[0] = 0;
  for (int i = 0; i < BlkSize; i++) {
    BlkStriIdx[i+1] = BlkStriIdx[i] + BlkStri[i];
  }
  free(BlkStri);

  /*
   * Allocate Memory - FaceVtx - FaceCell 
   */
  mesh->_dFaceVtx     = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t *) * dataSize/2); /* Not stupid at all */
  mesh->_dFaceVtxIdx  = (int *) malloc( sizeof(int *) * dataSize/2); /* Surdim as Hell */
  mesh->_dFaceCell    = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t *) * dataSize/2); /* Surdim as Hell */

  /*
   * Init AbsNumerotation 
   */
  int iAbsFace = 0;
  mesh->_dFaceVtxIdx[0] = 0;

  /*
   * Loop over BlkData : Each block correspond to a key 
   *              Each key can identify the same faces 
   *              Need to solve all conflict -> See _find_pairs
   */
  for(int iBlk=0; iBlk < BlkSize; iBlk++){

    /*
     * Traitement des groupes de données
     * Reminder Blk = [ [iCell, iType, nVtx, i1, i2, ...], ... ] 
     */
    int nEntryFace = BlkStriIdx[iBlk+1] - BlkStriIdx[iBlk];

    /* 
     * Verbose 
     */
    if(0 == 1){
      printf("BlkStri[%d] : %d -> %d \n", iBlk, BlkStriIdx[iBlk], nEntryFace);
      for(int j=BlkStriIdx[iBlk]; j < BlkStriIdx[iBlk+1]; j++){
        printf("BlkData[%d] : %d \n", j, BlkData[j]);
      }
    }

    /*
     * First Pass - Identify face conflict 
     */

    int nFac = 0;
    int iFacData = BlkStriIdx[iBlk];
    while(iFacData < BlkStriIdx[iBlk+1]){
      // printf("Number of Vextex for Face : %d \n", BlkData[iFacData+2]);
      iFacData += 3+BlkData[iFacData+2];
      nFac += 1;
    }

    /*
     *  Allocate array of displacement for all faces in conflict
     *         -> Set the adress of the block data of current faces
     *         -> Make a surdim before and not realloc ...
     */
    int *IdxFace = (int *) malloc( sizeof(int *) * nFac + 1);
    
    /** New **/
    iFacData   = BlkStriIdx[iBlk];
    IdxFace[0] = BlkStriIdx[iBlk];
    nFac       = 0;
    
    /** All Faces **/
    while(iFacData < BlkStriIdx[iBlk+1]){
      IdxFace[nFac+1] = IdxFace[nFac] + BlkData[iFacData+2] + 3;
      // printf("BlkData : %d \n ", BlkData[iFacData+2]);
      // printf("IdxFace[%d] : %d \n ", nFac+1, IdxFace[nFac+1]);
      iFacData     += BlkData[iFacData+2] + 3;
      nFac         += 1;
    }
    /* ************************************************** */
    

    /* 
     * Verbose 
     */
    if(0 == 1){
      for(int j=0; j < nFac+1; j++){
        printf("IdxFace[%d] : %d \n", j, IdxFace[j]);
        // printf("DeltaIdxFace[%d] : %d \n", j, IdxFace[j+1]-IdxFace[j]);
      }
    }

    /*
     * Solve conflict between all faces and build connectivity array 
     *   WATCH OUT -> iAbsFace is changed 
     */
    iAbsFace = _find_pairs(IdxFace, 
                           BlkData, 
                           nFac, 
                           iAbsFace, 
                           mesh->_dFaceVtx,
                           mesh->_dFaceVtxIdx,
                           mesh->_dFaceCell);

    /*
     * Free 
     */
    free(IdxFace);

  }  
  
  /* 
   * Fill up fetch structure
   */
  mesh->dNFace = iAbsFace;

  /* 
   * Realloc -> TODOUX
   */
  mesh->_dFaceVtxIdx = (int *        ) realloc((mesh->_dFaceVtxIdx), (mesh->dNFace + 1) * sizeof(int * ) );
  mesh->_dFaceVtx    = (PDM_g_num_t *) realloc((mesh->_dFaceVtx   ), mesh->_dFaceVtxIdx[mesh->dNFace] * sizeof(PDM_g_num_t * ));
  mesh->_dFaceCell   = (PDM_g_num_t *) realloc((mesh->_dFaceCell  ), mesh->dNFace * 2 * sizeof(PDM_g_num_t * ));

  /* 
   * Verbose
   */
  if(1 == 0){
    printf("fetch->dNFace : %d \n", mesh->dNFace);
    for(int i = 0; i < mesh->dNFace; i++) {
      // printf("BlockData[%d]    : %d \n", i, BlkData2[i]);
      printf("mesh->_dFaceCell[%d] : %d/%d \n", i, mesh->_dFaceCell[2*i], mesh->_dFaceCell[2*i+1]);
    }

  }

  /* 
   * Generate absolute numerotation of faces 
   */
  _make_absolute_face_numbering(mesh);
    
  /* 
   * Rebuild cell face 
   */
  PDM_g_num_t *LNToGNElem = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * 2 * mesh->dNFace );
  for (PDM_g_num_t i = mesh->face_distrib[mesh->i_proc]; i < mesh->face_distrib[mesh->i_proc+1]; i++) {
    int idx = (int) (i-mesh->face_distrib[mesh->i_proc]);
    LNToGNElem[2*idx  ] = i+1;
    LNToGNElem[2*idx+1] = i+1;
  }
  
  /* 
   * Verbose 
   */
  if(0 == 1){
    printf("------------------------------------------- \n");
    for (int i = mesh->face_distrib[mesh->i_proc]; i < mesh->face_distrib[mesh->i_proc+1]; i++) {
      // printf("[%d] - i/iFace  : %d - %d \n", mesh->i_proc, i, LNToGNElem[i]);
      printf("[%d] - i/iFace  : %d - %d / %d \n", mesh->i_proc, i, LNToGNElem[2*i], LNToGNElem[2*i+1]);
    }
    printf("------------------------------------------- \n");
  }
  
    
  int nFac2 = 2*mesh->dNFace;
  
  int* part_stri2 = (int *) malloc( sizeof(int) * 2 * mesh->dNFace );
  
  /* Prepare dFaceCell loc */
  PDM_g_num_t* dFaceCellTmp = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * 2 * mesh->dNFace );
  
  int idxG = 0;
  for (int i = mesh->face_distrib[mesh->i_proc]; i < mesh->face_distrib[mesh->i_proc+1]; i++) {
    int idx = (int) (i-mesh->face_distrib[mesh->i_proc]);
    dFaceCellTmp[idxG] = mesh->_dFaceCell[2*idx];
    LNToGNElem[idxG] = i+1;
    part_stri2[idxG] = 1;
    idxG++;
    if(mesh->_dFaceCell[2*idx+1] != 0){
      dFaceCellTmp[idxG] = mesh->_dFaceCell[2*idx+1];
      LNToGNElem[idxG  ] = i+1;
      part_stri2[idxG] = 1;
      idxG++;
    }
    else
    {
      dFaceCellTmp[idxG] = mesh->_dFaceCell[2*idx];
      LNToGNElem[idxG  ] = i+1;
      part_stri2[idxG] = 0;
      idxG++;
    }
  }
  
  nFac2 = idxG;
  
  
  
  PDM_part_to_block_t *ptb2 = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1.,
                                                       &dFaceCellTmp,
                                                       // &mesh->_dFaceCell,
                                                       &nFac2,
                                                       1,
                                                       mesh->pdm_mpi_comm); 

  int *BlkStri2 = NULL;
  int *BlkData2 = NULL;
  
  /** Prepare part_stride for dFaceCell - 2 for all **/
  // for(int i = 0; i < 2*mesh->dNFace ; i++){
  // for(int i = 0; i < idxG ; i++){
  //   part_stri2[i] = 1;
  //   // part_stri2[i] = mesh->dNFace;
  // }
    
  PDM_part_to_block_exch(          ptb2,
                                   sizeof(PDM_g_num_t),
                                   PDM_STRIDE_VAR,
                                   -1,
                                   &part_stri2,
                         (void **) &LNToGNElem,
                                   &BlkStri2,
                         (void **) &BlkData2);

  /*
   *  Get the size of the current process bloc
   */
  int dElmtTot = PDM_part_to_block_n_elt_block_get(ptb2); /* Volumic and surfacic */
  mesh->n_dcell = dElmtTot;

  /* 
   * Verbose 
   */
  if(0 == 1){
    printf("dElmtTot : %d\n", dElmtTot);
    for(int i = 0; i < dElmtTot; i++) {
      printf("BlockData[%d]    : %d \n", i, BlkData2[i]);
      printf("BlkStri2[%d] : %d \n", i, BlkStri2[i]);
    }
  }
  
  mesh->dCellFaceIdx = (PDM_l_num_t * ) malloc( (dElmtTot + 1) * sizeof(PDM_l_num_t * ) );
  
  mesh->dCellFaceIdx[0] = 0;
  for(int i = 0; i < dElmtTot; i++){
    mesh->dCellFaceIdx[i+1] = mesh->dCellFaceIdx[i] + BlkStri2[i];
  }
  
  // printf("dElmtTot : %d\n", mesh->dCellFaceIdx[dElmtTot]);
  
  mesh->dCellFace = (PDM_g_num_t *) malloc(mesh->dCellFaceIdx[dElmtTot] * sizeof(PDM_g_num_t * ));
   
  for(int i = 0; i < dElmtTot; i++){
    
    // int nFacPerElmt = BlkStri2[i];
    int nFacPerElmt = mesh->dCellFaceIdx[i+1] - mesh->dCellFaceIdx[i];
    // printf("nFacPerElmt : %d\n", nFacPerElmt);
    for(int iFac = 0; iFac < nFacPerElmt; iFac++){
      mesh->dCellFace[mesh->dCellFaceIdx[i]+iFac] = BlkData2[mesh->dCellFaceIdx[i]+iFac];
    }
  }
  
  /* 
   * Mask all FaceCell with boundary 
   */
  // for(int i = 0; i < mesh->dNFace; i++) {
  //    if(mesh->_dFaceCell[2*i+1] > mesh->n_cell_abs){
  //       mesh->_dFaceCell[2*i+1] = 0;
  //    }
  // }
  
  /* Free */
  PDM_part_to_block_free(ptb );
  PDM_part_to_block_free(ptb2);

  free(LNToGN);
  free(LNToGNElem);
  free(dFaceCellTmp);
  free(part_data);
  free(part_stri);

  free(BlkData   );
  free(BlkStriIdx);
  free(BlkStri2  );
  free(BlkData2  );
  
  // PDM_error (__FILE__, __LINE__, 0, "Not implemented yet\n");
}


/**
 * \brief  Return cell \rightarrow face connectivity
 *
 * \param [in]   hdl              Distributed nodal mesh handle
 * \param [out]  dCellFaceIdx   Index of distributed cell->face connectivity
 * \param [out]  dCellFace       Distributed cell->face connectivity
 *
 * \return     Number of cells on the current process
 *  
 */

int
PDM_DMesh_nodal_cell_face_get
(
const int   hdl,
      int   **dCellFaceIdx,  
PDM_g_num_t **dCellFace  
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  *dCellFaceIdx = mesh->dCellFaceIdx;  
  *dCellFace = mesh->dCellFace;

  return mesh->n_dcell;

}

/**
 * \brief  Return face->cell connectivity
 *
 * \param [in]   hdl              Distributed nodal mesh handle
 * \param [out]  FaceCell       Distributed face->cell connectivity
 *
 * \return     Number of cells on the current process
 *  
 */
int
PDM_DMesh_nodal_face_cell_get
(
const int     hdl,
PDM_g_num_t **dFaceCell  
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  *dFaceCell = mesh->_dFaceCell;

  return mesh->dNFace;
  
}



/**
 * \brief  Return face \rightarrow vertex connectivity
 *
 * \param [in]   hdl              Distributed nodal mesh handle
 * \param [out]  dCellFaceIdx   Index of distributed cell->face connectivity
 * \param [out]  dCellFace       Distributed cell->face connectivity
 *
 * \return     Number of faces on the current process
 *  
 */

int
PDM_DMesh_nodal_face_vtx_get
(
const int   hdl,
      int   **_dFaceVtxIdx,  
PDM_g_num_t **_dFaceVtx  
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  *_dFaceVtxIdx = mesh->_dFaceVtxIdx;  
  *_dFaceVtx     = mesh->_dFaceVtx;

  return mesh->dNFace; 

}

/**
 * \brief  Return vertices distribution
 * 
 * \param [in]  hdl  Distributed nodal mesh handle
 *
 * \return  A array of size \ref n_procs + 1    
 * 
 */

const PDM_g_num_t *
PDM_DMesh_nodal_distrib_vtx_get
(
 const int          hdl
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }
  
  PDM_DMesh_nodal_vtx_t *vtx = mesh->vtx;

  return vtx->distrib;

}


/**
 * \brief  Return section distribution
 * 
 * \param [in]  hdl        Distributed nodal mesh handle
 * \param [in]  id_section   Block identifier
 *
 * \return  A array of size \ref n_procs + 1    
 * 
 */

const PDM_g_num_t *
PDM_DMesh_nodal_distrib_section_get
(
 const int   hdl,
 const int   id_section     
)
{
  PDM_DMesh_nodal_t *mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  int _id_section;
  
  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {
  
    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;
  
    const PDM_DMesh_nodal_section_poly3d_t *section = 
    (const PDM_DMesh_nodal_section_poly3d_t *) 
     PDM_Handles_get (mesh->sections_poly3d, _id_section);
  
    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
    }
      
    return section->distrib;
  }
  
  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {
  
    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
  
    const PDM_DMesh_nodal_section_poly2d_t *section = 
    (const PDM_DMesh_nodal_section_poly2d_t *) 
     PDM_Handles_get (mesh->sections_poly2d, _id_section);
  
    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad polygon section identifier\n");
    }
  
    return section->distrib;
  }
  
  else {
  
    _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
  
    const PDM_DMesh_nodal_section_std_t *section = 
    (const PDM_DMesh_nodal_section_std_t *) 
     PDM_Handles_get (mesh->sections_std, _id_section);
  
    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad polyhedron section identifier\n");
    }
  
    return section->distrib;
  }
}


/**
 * \brief  Return cell distribution
 * 
 * \param [in]  hdl  Distributed nodal mesh handle
 *
 * \return  A array of size \ref n_procs + 1    
 * 
 */

const PDM_g_num_t *
PDM_DMesh_nodal_distrib_cell_get
(
 const int  hdl
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  return mesh->cell_distrib;
  
}


/**
 * \brief  Return face distribution
 * 
 * \param [in]  hdl  Distributed nodal mesh handle
 *
 * \return  A array of size \ref n_procs + 1    
 * 
 */

const PDM_g_num_t *
PDM_DMesh_nodal_distrib_face_get
(
 const int hdl
)
{
  PDM_DMesh_nodal_t * mesh = 
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);
  
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");  
  }

  return mesh->face_distrib;
  
}


#ifdef __cplusplus
}
#endif /* __cplusplus */

