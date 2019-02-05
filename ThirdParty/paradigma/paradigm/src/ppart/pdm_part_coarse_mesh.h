/* 
 * File:   pdm_part_coarse_mesh.h
 * Author: jmagnene
 *
 * Created on July 8, 2016, 9:29 AM
 */

#ifndef PDM_PART_COARSE_MESH_H
#define	PDM_PART_COARSE_MESH_H


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_part.h"


#if !defined (__hpux) && !defined (_AIX) 
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif


#ifdef	__cplusplus
extern "C" {
#endif


/**
 *
 * \brief Return an initialized coarse mesh object
 *
 * \param [out]  cmId              Coarse mesh identifier
 * 
 * \param [in]   pt_comm           Communicator
 * \param [in]   method            Choice between (1 for ParMETIS or 2 for PT-Scotch)
 * \param [in]   nPart             Number of partitions
 * \param [in]   nTPart            Total number of partitions
 * \param [in]   have_cellTag      Presence d'un tableau de tags pour les cellules
 * \param [in]   have_faceTag      Presence d'un tableau de tags pour les faces
 * \param [in]   have_vtxTag       Presence d'un tableau de tags pour les sommets
 * \param [in]   have_cellWeight   Presence d'un tableau de poids pour les cellules
 * \param [in]   have_faceWeight   Presence d'un tableau de poids pour les faces
 * \param [in]   have_faceGroup    Presence des tableaux de groupes de faces
 */
    
void 
PDM_part_coarse_mesh_create
(
 int                *cmId,
 PDM_MPI_Comm        comm,        
 const char*         method,
 const int           nPart, 
 const int           nTPart,
 const int           nFaceGroup,
 const int           have_cellTag,
 const int           have_faceTag,
 const int           have_vtxTag,
 const int           have_cellWeight,
 const int           have_faceWeight,
 const int           have_faceGroup
);

void
PROCF (pdm_part_coarse_mesh_create_cf, PDM_PART_COARSE_MESH_CREATE_CF)
(
 int                *cmId,
 PDM_MPI_Fint       *fcomm,        
 const char         *method,
 const int          *l_method,
 const int          *nPart, 
 const int          *nTPart, 
 const int          *nFaceGroup,
 const int          *have_cellTag,
 const int          *have_faceTag,
 const int          *have_vtxTag,
 const int          *have_cellWeight, 
 const int          *have_faceWeight,
 const int          *have_faceGroup
 );


/**
 *
 * \brief Build a coarse mesh
 *
 * \param [in]  cmId               Coarse mesh identifier 
 * \param [in]  iPart              Partition identifier
 * \param [in]  nCoarseCell        Number of cells in the coarse grid
 * \param [in]  nCell              Number of cells
 * \param [in]  nFace              Number of faces
 * \param [in]  nFacePartBound     Number of partitioning boundary faces
 * \param [in]  nVtx               Number of vertices 
 * \param [in]  nFaceGroup         Number of face groups             
 * \param [in]  cellFaceIdx        Cell to face connectivity index (size = nCell + 1, numbering : 0 to n-1)
 * \param [in]  cellFace           Cell to face connectivity (size = cellFaceIdx[nCell] = lCellFace
 *                                                             numbering : 1 to n)
 * \param [in]  cellTag            Cell tag (size = nCell)
 * \param [in]  cellLNToGN         Cell local numbering to global numbering (size = nCell, numbering : 1 to n)
 * \param [in]  cellWeight         Cell weight (size = nCell)
 * \param [in]  faceWeight         Face weight (size = nFace)
 * \param [in]  faceCell           Face to cell connectivity  (size = 2 * nFace, numbering : 1 to n)
 * \param [in]  faceVtxIdx         Face to Vertex connectivity index (size = nFace + 1, numbering : 0 to n-1)
 * \param [in]  faceVtx            Face to Vertex connectivity (size = faceVertexIdx[nFace], numbering : 1 to n)
 * \param [in]  faceTag            Face tag (size = nFace)
 * \param [in]  faceLNToGN         Face local numbering to global numbering (size = nFace, numbering : 1 to n)
 * \param [in]  vtxCoord           Vertex coordinates (size = 3 * nVertex)
 * \param [in]  vtxTag             Vertex tag (size = nVertex)
 * \param [in]  vtxLNToGN          Vertex local numbering to global numbering (size = nVtx, numbering : 1 to n)
 * \param [in]  faceGroupIdx       Face group index (size = nFaceGroup + 1, numbering : 1 to n-1)
 * \param [in]  faceGroup          faces for each group (size = faceGroupIdx[nFaceGroup] = lFaceGroup, numbering : 1 to n)
 * \param [in]  faceGroupLNToGN    Faces global numbering for each group 
 *                                  (size = faceGroupIdx[nFaceGroup] = lFaceGroup, numbering : 1 to n)
 * \param [in]  facePartBoundProcIdx  Partitioning boundary faces block distribution from processus (size = nProc + 1)
 * \param [in]  facePartBoundPartIdx  Partitioning boundary faces block distribution from partition (size = nTPart + 1)
 * \param [in]  facePartBound      Partitioning boundary faces (size = 4 * nFacePartBound)
 *                                       sorted by processus, sorted by partition in each processus, and
 *                                       sorted by absolute face number in each partition
 *                                   For each face :
 *                                        - Face local number (numbering : 1 to n)
 *                                        - Connected process (numbering : 0 to n-1)
 *                                        - Connected Partition 
 *                                          on the connected process (numbering :1 to n)
 *                                        - Connected face local number 
 *                                          in the connected partition (numbering :1 to n)
 */

void 
PDM_part_coarse_mesh_input
(
 int                 cmId,
 int                 iPart,
 const int           nCoarseCell,
 const int           nCell,
 const int           nFace,
 const int           nVtx,
 const int           nFaceGroup,
 const int           nFacePartBound,
 const int          *cellFaceIdx,
 const int          *cellFace,
 const int          *cellTag,
 const int          *cellWeight,
 const int          *faceWeight,
 const PDM_g_num_t  *cellLNToGN,       
 const int          *faceCell,
 const int          *faceVtxIdx,
 const int          *faceVtx,
 const int          *faceTag,       
 const PDM_g_num_t  *faceLNToGN,       
 const double       *vtxCoord,
 const int          *vtxTag,
 const PDM_g_num_t  *vtxLNToGN,       
 const int          *faceGroupIdx,
 const int          *faceGroup,
 const PDM_g_num_t  *faceGroupLNToGN,
 const int          *facePartBoundProcIdx,       
 const int          *facePartBoundPartIdx,
 const int          *facePartBound
);

void
PROCF (pdm_part_coarse_mesh_input, PDM_PART_COARSE_MESH_INPUT)
(
 int                *cmId,
 int                *iPart,
 const int          *nCoarseCell,
 const int          *nCell,
 const int          *nFace,
 const int          *nVtx,
 const int          *nFaceGroup,
 const int          *nFacePartBound,
 const int          *cellFaceIdx,
 const int          *cellFace,
 const int          *have_cellTag,
 const int          *cellTag,
 const int          *have_cellWeight,
 const int          *cellWeight,
 const int          *have_faceWeight,
 const int          *faceWeight,
 const PDM_g_num_t *cellLNToGN,       
 const int          *faceCell,
 const int          *faceVtxIdx,
 const int          *faceVtx,
 const int          *have_faceTag,
 const int          *faceTag,       
 const PDM_g_num_t *faceLNToGN,       
 const double       *vtxCoord,
 const int          *have_vtxTag,
 const int          *vtxTag,
 const PDM_g_num_t *vtxLNToGN,       
 const int          *have_faceGroup,
 const int          *faceGroupIdx,
 const int          *faceGroup,
 const PDM_g_num_t *faceGroupLNToGN,
 const int          *facePartBoundProcIdx,       
 const int          *facePartBoundPartIdx,
 const int          *facePartBound
);
 
/**
 *
 * \brief Updates all the arrays dealing with MPI exchanges
 *
 * \param [in] cmId               Coarse mesh identifier
 */

void 
PDM_part_coarse_mesh_compute
(
  int cmId
);

void
PROCF (pdm_part_coarse_mesh_compute, PDM_PART_COARSE_MESH_COMPUTE)
(
  int *cmId
);

/**
 *
 * \brief Return a coarse mesh partition dimensions
 * 
 * \param [in]   cmId               Coarse mesh identifier
 * \param [in]   iPart              Current partition
 * 
 * \param [out]  nCell              Number of cells
 * \param [out]  nFace              Number of faces
 * \param [out]  nFacePartBound     Number of partitioning boundary faces
 * \param [out]  nVtx               Number of vertices 
 * \param [out]  nProc              Number of processus
 * \param [out]  nTPart             Number of partitions
 * \param [out]  sCellFace          Size of cell-face connectivity 
 * \param [out]  sFaceVtx           Size of face-vertex connectivity
 * \param [out]  sFacePartBound     Size of facePartBound array
 * \param [out]  sFaceGroup         Size of faceGroup array  
 * \param [out]  sCoarseCellToFineCell  Size of coarseCellToFineCell array
 *
 */

void 
PDM_part_coarse_mesh_part_dim_get
(
 const int      cmId,
 int            iPart,
 int           *nCell,
 int           *nFace,
 int           *nFacePartBound,
 int           *nVtx,
 int           *nProc,
 int           *nTPart,
 int           *nFaceGroup,
 int           *sCellFace,
 int           *sFaceVtx,
 int           *sFaceGroup,
 int           *sCoarseCellToFineCell
);

void
PROCF (pdm_part_coarse_mesh_part_dim_get, PDM_PART_COARSE_MESH_PART_DIM_GET)
(
 const int     *cmId,
 int           *iPart,
 int           *nCell,
 int           *nFace,
 int           *nFacePartBound,
 int           *nVtx,
 int           *nProc,
 int           *nTPart,
 int           *nFaceGroup,
 int           *sCellFace,
 int           *sFaceVtx,
 int           *sFaceGroup,
 int           *sCoarseCellToFineCell
);

/**
 *
 * \brief Return a mesh partition
 * 
 * \param [in]   cmId               Coarse mesh identifier
 * \param [in]   iPart              Current partition
 * 
 * \param [out]  cellFaceIdx        Cell to face connectivity index (size = nCell + 1, numbering : 0 to n-1)
 * \param [out]  cellFace           Cell to face connectivity (size = cellFaceIdx[nCell] = lCellFace
 *                                                             numbering : 1 to n)
 * \param [out]  cellTag            Cell tag (size = nCell)
 * \param [out]  cellLNToGN         Cell local numbering to global numbering (size = nCell, numbering : 1 to n)
 * \param [out]  cellInitCellIdx    Array of indexes of the connected partitions (size : nCoarseCell + 1)
 * \param [out]  cellInitCell       Partitioning array (size : cellInitCellIdx[nCoarseCell]) 
 * 
 * \param [out]  faceCell           Face to cell connectivity  (size = 2 * nFace, numbering : 1 to n)
 * \param [out]  faceVtxIdx         Face to Vertex connectivity index (size = nFace + 1, numbering : 0 to n-1)
 * \param [out]  faceVtx            Face to Vertex connectivity (size = faceVertexIdx[nFace], numbering : 1 to n)
 * \param [out]  faceTag            Face tag (size = nFace)
 * \param [out]  faceLNToGN         Face local numbering to global numbering (size = nFace, numbering : 1 to n)
 * \param [out]  faceInitFace       Coarse face - fine face connectivity (size = nCoarseFace)
 * 
 * \param [out]  vtxCoord           Vertex coordinates (size = 3 * nVtx)
 * \param [out]  vtxTag             Vertex tag (size = nVtx)
 * \param [out]  vtxLNToGN          Vertex local numbering to global numbering (size = nVtx, numbering : 1 to n)
 * \param [out]  vtxInitVtx         Coarse vertex - fine vertex connectivity (size = nCoarseVtx)
 * 
 * \param [out]  faceGroupIdx       Face group index (size = nFaceGroup + 1, numbering : 1 to n-1)
 * \param [out]  faceGroup          faces for each group (size = faceGroupIdx[nFaceGroup] = lFaceGroup, numbering : 1 to n)
 * \param [out]  faceGroupLNToGN    Faces global numbering for each group 
 *                                  (size = faceGroupIdx[nFaceGroup] = lFaceGroup, numbering : 1 to n)
 * 
 * \param [out]  facePartBoundProcIdx  Partitioning boundary faces block distribution from processus (size = nProc + 1)
 * \param [out]  facePartBoundPartIdx  Partitioning boundary faces block distribution from partition (size = nTPart + 1)
 * \param [out]  facePartBound      Partitioning boundary faces (size = 4 * nFacePartBound)
 *                                       sorted by processus, sorted by partition in each processus, and
 *                                       sorted by absolute face number in each partition
 *                                   For each face :
 *                                        - Face local number (numbering : 1 to n)
 *                                        - Connected process (numbering : 0 to n-1)
 *                                        - Connected Partition 
 *                                          on the connected process (numbering :1 to n)
 *                                        - Connected face local number 
 *                                          in the connected partition (numbering :1 to n)
 * 
 */

void 
PDM_part_coarse_mesh_part_get
(
 const int    cmId,
 const int    iPart,       
 int          **cellFaceIdx,
 int          **cellFace,
 int          **cellTag,
 PDM_g_num_t **cellLNToGN,
 int          **cellInitCellIdx,                  
 int          **cellInitCell,          
 int          **faceCell,
 int          **faceVtxIdx,
 int          **faceVtx,
 int          **faceTag,
 PDM_g_num_t **faceLNToGN,  
 int          **faceGroupInitFaceGroup,
 int          **faceInitFace,          
 double       **vtxCoord,
 int          **vtxTag,
 PDM_g_num_t **vtxLNToGN,        
 int          **vtxInitVtx,          
 int          **faceGroupIdx,
 int          **faceGroup,
 PDM_g_num_t **faceGroupLNToGN,
 int          **facePartBoundProcIdx,
 int          **facePartBoundPartIdx,
 int          **facePartBound        
);

void
PROCF (pdm_part_coarse_mesh_part_get, PDM_PART_COARSE_MESH_PART_GET)
(
 int          *cmId,
 int          *iPart,       
 int          *cellFaceIdx,
 int          *cellFace,
 int          *cellTag,
 PDM_g_num_t  *cellLNToGN,
 int          *cellInitCellIdx,                  
 int          *cellInitCell,          
 int          *faceCell,
 int          *faceVtxIdx,
 int          *faceVtx,
 int          *faceTag,
 PDM_g_num_t *faceLNToGN,  
 int          *faceGroupInitFaceGroup,
 int          *faceInitFace,          
 double       *vtxCoord,
 int          *vtxTag,
 PDM_g_num_t *vtxLNToGN,        
 int          *vtxInitVtx,          
 int          *faceGroupIdx,
 int          *faceGroup,
 PDM_g_num_t *faceGroupLNToGN,
 int          *facePartBoundProcIdx,
 int          *facePartBoundPartIdx,
 int          *facePartBound
);

/**
 *
 * \brief Free coarse mesh
 *
 * \param [in]   cmId        Coarse mesh identifier
 *
 */

void 
PDM_part_coarse_mesh_free
(
 const int    cmId
);

void
PROCF (pdm_part_coarse_mesh_free, PDM_PART_COARSE_MESH_FREE)
(
  int *cmId
);

/**
 *
 * \brief Return times
 * 
 * \param [in]   cmId        coarse mesh identifier
 * \param [out]  elapsed     elapsed times (size = 18)
 * \param [out]  cpu         cpu times (size = 18)
 * \param [out]  cpu_user    user cpu times (size = 18)
 * \param [out]  cpu_sys     system cpu times (size = 18)
 *
 */

void PDM_part_coarse_mesh_time_get
(
 int       cmId,
 double  **elapsed,
 double  **cpu,
 double  **cpu_user,
 double  **cpu_sys
);

void 
PROCF (pdm_part_coarse_mesh_time_get, PDM_PART_COARSE_MESH_TIME_GET)
(
 int      *cmId,
 double   *elapsed,
 double   *cpu,
 double   *cpu_user,
 double   *cpu_sys
 );

/**
 *
 * \brief Displays all the arrays of a coarse mesh
 * 
 * \param [in]   cmId        Coarse mesh identifier
 * 
 */

void 
PDM_part_coarse_mesh_display
(
 const int    cmId
);

void  
PROCF (pdm_part_coarse_mesh_display, PDM_PART_COARSE_MESH_DISPLAY)
(
  int *cmId
);
  

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_PART_COARSE_MESH_H */

