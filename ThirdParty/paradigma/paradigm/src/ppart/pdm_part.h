#ifndef __PDM_PART_H__
#define __PDM_PART_H__

#include <stdio.h>
#include "pdm.h"
#include "pdm_mpi.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \enum PDM_part_split_t
 * \brief Split method
 *
 */

typedef enum {
  PDM_PART_SPLIT_PARMETIS = 1,
  PDM_PART_SPLIT_PTSCOTCH = 2,
  PDM_PART_SPLIT_HILBERT  = 3
} PDM_part_split_t;

typedef struct _PDM_part_t PDM_part_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Build a initial partitioning
 *
 *  Build a initial partitioning from :
 *      - Cell block distribution with implicit global numbering
 *         (the first cell is the first cell of the first process and
 *          the latest cell is the latest cell of the latest process)
 *      - Face block distribution with implicit global numbering
 *      - Vertex block distribution with implicit global numbering
 *  To repart an existing partition use \ref PDM_part_repart function
 *
 * \param [out]  ppartId        ppart identifier
 * \param [in]   pt_comm        MPI Comminicator
 * \param [in]   split_method   Split method
 * \param [in]   renum_cell_method Cell renumbering method
 * \param [in]   renum_face_method Cell renumbering method
 * \param [in]   renum_properties_cell  For cache blocking [ nCellPerCacheWanted, isAsynchrone, isVectorisation ] \ref PDM_renum_cacheblocking
 * \param [in]   renum_face_method Cell renumbering method
 * \param [in]   renum_properties_face  NOT USE
 * \param [in]   nPart          Number of partition to build on this process
 * \param [in]   dNCell         Number of distributed cells
 * \param [in]   dNFace         Number of distributed faces
 * \param [in]   dNVtx          Number of distributed vertices
 * \param [in]   nFaceGroup     Number of face groups
 * \param [in]   dCellFaceIdx   Distributed cell face connectivity index or NULL
 *                              (size : dNCell + 1, numbering : 0 to n-1)
 * \param [in]   dCellFace      Distributed cell face connectivity or NULL
 *                              (size : dFaceVtxIdx[dNCell], numbering : 1 to n)
 * \param [in]   dCellTag       Cell tag (size : nCell) or NULL
 * \param [in]   dCellWeight    Cell weight (size : nCell) or NULL
 * \param [in]   dCellPart      Distributed cell partitioning
 *                              (size = dNCell) or NULL (No partitioning if != NULL)
 * \param [in]   dFaceCell      Distributed face cell connectivity or NULL
 *                              (size : 2 * dNFace, numbering : 1 to n)
 * \param [in]   dFaceVtxIdx    Distributed face to vertex connectivity index
 *                              (size : dNFace + 1, numbering : 0 to n-1)
 * \param [in]   dFaceVtx       Distributed face to vertex connectivity
 *                              (size : dFaceVtxIdx[dNFace], numbering : 1 to n)
 * \param [in]   dFaceTag       Distributed face tag (size : dNFace)
 *                              or NULL
 * \param [in]   dVtxCoord      Distributed vertex coordinates
 *                              (size : 3*dNVtx)
 * \param [in]   dVtxTag        Distributed vertex tag (size : dNVtx) or NULL
 * \param [in]   dFaceGroupIdx  Index of distributed faces list of each group
 *                              (size = nFaceGroup + 1) or NULL
 * \param [in]   dFaceGroup     distributed faces list of each group
 *                              (size = dFaceGroup[dFaceGroupIdx[nFaceGroup]], numbering : 1 to n)
 *                              or NULL
 *
 */

void
PDM_part_create
(
 int                         *ppartId,
 const PDM_MPI_Comm           comm,
 const PDM_part_split_t       split_method,
 const char                  *renum_cell_method,
 const char                  *renum_face_method,
 const int                    nPropertyCell,
 const int                   *renum_properties_cell,
 const int                    nPropertyFace,
 const int                   *renum_properties_face,
 const int                    nPart,
 const int                    dNCell,
 const int                    dNFace,
 const int                    dNVtx,
 const int                    nFaceGroup,
 const int                   *dCellFaceIdx,
 const PDM_g_num_t           *dCellFace,
 const int                   *dCellTag,
 const int                   *dCellWeight,
 const int                    have_dCellPart,
       int                   *dCellPart,
 const PDM_g_num_t           *dFaceCell,
 const int                   *dFaceVtxIdx,
 const PDM_g_num_t           *dFaceVtx,
 const int                   *dFaceTag,
 const double                *dVtxCoord,
 const int                   *dVtxTag,
 const int                   *dFaceGroupIdx,
 const PDM_g_num_t           *dFaceGroup
 );

void
PROCF (pdm_part_create_cf, PDM_PART_CREATE_CF)
(
 int                *ppartId,
 const PDM_MPI_Fint *fcomm,
 const int          *split_method,
 const char         *renum_cell_method,
 const int          *l_renum_cell_method,
 const char         *renum_face_method,
 const int          *l_renum_face_method,
 const int          *nPropertyCell,
 const int          *renum_properties_cell,
 const int          *nPropertyFace,
 const int          *renum_properties_face,
 const int          *nPart,
 const int          *dNCell,
 const int          *dNFace,
 const int          *dNVtx,
 const int          *nFaceGroup,
 const int          *have_dCellFace,
 const int          *dCellFaceIdx,
 const PDM_g_num_t  *dCellFace,
 const int          *have_dCellTag,
 const int          *dCellTag,
 const int          *have_dCellWeight,
 const int          *dCellWeight,
 const int          *have_dCellPart,
       int          *dCellPart,
 const int          *have_dFaceCell,
 const PDM_g_num_t  *dFaceCell,
 const int          *dFaceVtxIdx,
 const PDM_g_num_t  *dFaceVtx,
 const int          *have_dFaceTag,
 const int          *dFaceTag,
 const double       *dVtxCoord,
 const int          *have_dVtxTag,
 const int          *dVtxTag,
 const int          *dFaceGroupIdx,
 const PDM_g_num_t  *dFaceGroup
);

/**
 *
 * \brief Return a mesh partition dimensions
 *
 * \param [in]   ppartId            ppart identifier
 * \param [in]   ipart              Current partition
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
 *
 */

void
PDM_part_part_dim_get
(
const  int    ppartId,
const  int    ipart,
       int   *nCell,
       int   *nFace,
       int   *nFacePartBound,
       int   *nVtx,
       int   *nProc,
       int   *nTPart,
       int   *sCellFace,
       int   *sFaceVtx,
       int   *sFaceGroup,
       int   *nFaceGroup
);

void
PROCF (pdm_part_part_dim_get, PDM_PART_PART_DIM_GET)
(
 int           *ppartId,
 int           *ipart,
 int           *nCell,
 int           *nFace,
 int           *nFacePartBound,
 int           *nVtx,
 int           *nProc,
 int           *nTPart,
 int           *sCellFace,
 int           *sFaceVtx,
 int           *sFaceGroup,
 int           *nFaceGroup
);

/**
 *
 * \brief Return a mesh partition
 *
 * \param [in]   ppartId               ppart identifier
 * \param [in]   ipart                 Current partition
 * \param [out]  cellTag               Cell tag (size = nCell)
 * \param [out]  cellFaceIdx           Cell to face connectivity index (size = nCell + 1, numbering : 0 to n-1)
 * \param [out]  cellFace              Cell to face connectivity (size = cellFaceIdx[nCell] = lCellFace
 *                                                                numbering : 1 to n)
 * \param [out]  cellLNToGN            Cell local numbering to global numbering (size = nCell, numbering : 1 to n)
 * \param [out]  faceTag               Face tag (size = nFace)
 * \param [out]  faceCell              Face to cell connectivity  (size = 2 * nFace, numbering : 1 to n)
 * \param [out]  faceVtxIdx            Face to Vertex connectivity index (size = nFace + 1, numbering : 0 to n-1)
 * \param [out]  faceVtx               Face to Vertex connectivity (size = faceVertexIdx[nFace], numbering : 1 to n)
 * \param [out]  faceLNToGN            Face local numbering to global numbering (size = nFace, numbering : 1 to n)
 * \param [out]  facePartBoundProcIdx  Partitioning boundary faces block distribution from processus (size = nProc + 1)
 * \param [out]  facePartBoundPartIdx  Partitioning boundary faces block distribution from partition (size = nTPart + 1)
 * \param [out]  facePartBound         Partitioning boundary faces (size = 4 * nFacePartBound)
 *                                          sorted by processus, sorted by partition in each processus, and
 *                                          sorted by absolute face number in each partition
 *                                      For each face :
 *                                           - Face local number (numbering : 1 to n)
 *                                           - Connected process (numbering : 0 to n-1)
 *                                           - Connected Partition
 *                                             on the connected process (numbering :1 to n)
 *                                           - Connected face local number
 *                                             in the connected partition (numbering :1 to n)
 * \param [out]  vtxTag                Vertex tag (size = nVertex)
 * \param [out]  vtx                   Vertex coordinates (size = 3 * nVertex)
 * \param [out]  vtxLNToGN             Vertex local numbering to global numbering (size = nVtx, numbering : 1 to n)
 * \param [out]  faceGroupIdx          Face group index (size = nFaceGroup + 1, numbering : 1 to n-1)
 * \param [out]  faceGroup             faces for each group (size = faceGroupIdx[nFaceGroup] = lFaceGroup, numbering : 1 to n)
 * \param [out]  faceGroupLNToGN       Faces global numbering for each group
 *                                     (size = faceGroupIdx[nFaceGroup] = lFaceGroup, numbering : 1 to n)
 *
 */

void PDM_part_part_val_get
(
const int            ppartId,
const int            ipart,
      int          **cellTag,
      int          **cellFaceIdx,
      int          **cellFace,
      PDM_g_num_t  **cellLNToGN,
      int          **faceTag,
      int          **faceCell,
      int          **faceVtxIdx,
      int          **faceVtx,
      PDM_g_num_t  **faceLNToGN,
      int          **facePartBoundProcIdx,
      int          **facePartBoundPartIdx,
      int          **facePartBound,
      int          **vtxTag,
      double       **vtx,
      PDM_g_num_t  **vtxLNToGN,
      int          **faceGroupIdx,
      int          **faceGroup,
      PDM_g_num_t  **faceGroupLNToGN
);

void
PROCF (pdm_part_part_val_get, PDM_PART_PART_VAL_GET)
(
 int           *ppartId,
 int           *ipart,
 int           *cellTag,
 int           *cellFaceIdx,
 int           *cellFace,
 PDM_g_num_t   *cellLNToGN,
 int           *faceTag,
 int           *faceCell,
 int           *faceVtxIdx,
 int           *faceVtx,
 PDM_g_num_t   *faceLNToGN,
 int           *facePartBoundProcIdx,
 int           *facePartBoundPartIdx,
 int           *facePartBound,
 int           *vtxTag,
 double        *vtx,
 PDM_g_num_t   *vtxLNToGN,
 int           *faceGroupIdx,
 int           *faceGroup,
 PDM_g_num_t   *faceGroupLNToGN
);

/**
 *
 * \brief Return a mesh partition
 *
 * \param [in]   ppartId               ppart identifier
 * \param [in]   ipart                 Current partition
 * \param [out]  cellColor             Cell Color (size = nCell)
 * \param [out]  faceColor             Face Color (size = nFace)
 */

void PDM_part_part_color_get
(
const int            ppartId,
const int            ipart,
      int          **cellColor,
      int          **faceColor,
      int          **threadColor,
      int          **hyperPlaneColor
);

void
PROCF (pdm_part_part_color_get, PDM_PART_PART_COLOR_GET)
(
 int           *ppartId,
 int           *ipart,
 int           *cellColor,
 int           *faceColor,
 int           *threadColor,
 int           *hyperPlaneColor
);

/**
 *
 * \brief Free ppart
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

void
PDM_part_free
(
const  int                ppartId
);

void
PROCF (pdm_part_free, PDM_PART_FREE)
(
 int                *ppartId
);

/**
 *
 * \brief Return times
 *
 * \param [in]   ppartId     ppart identifier
 * \param [out]  elapsed     elapsed times (size = 4)
 * \param [out]  cpu         cpu times (size = 4)
 * \param [out]  cpu_user    user cpu times (size = 4)
 * \param [out]  cpu_sys     system cpu times (size = 4)
 *
 */

void
PDM_part_time_get
(
const int       ppartId,
      double  **elapsed,
      double  **cpu,
      double  **cpu_user,
      double  **cpu_sys
);


void
PROCF (pdm_part_time_get, PDM_PART_TIME_GET)
(
 int      *ppartId,
 double   *elapsed,
 double   *cpu,
 double   *cpu_user,
 double   *cpu_sys
);


/**
 *
 * \brief Return statistic
 *
 * \param [in]   ppartId                        ppart identifier
 * \param [out]  cells_average                  average of cells number
 * \param [out]  cells_median                   median of cells number
 * \param [out]  cells_std_deviation            standard deviation of cells number
 * \param [out]  cells_min                      minimum of cells nummber
 * \param [out]  cells_max                      maximum of cells nummber
 * \param [out]  bound_part_faces_average       average of partitioning boundary faces
 * \param [out]  bound_part_faces_median        median of partitioning boundary faces
 * \param [out]  bound_part_faces_std_deviation standard deviation of partitioning boundary faces
 * \param [out]  bound_part_faces_min           minimum of partitioning boundary faces
 * \param [out]  bound_part_faces_max           maximum of partitioning boundary faces
 *
 */

void
PDM_part_stat_get
(
const int       ppartId,
      int      *cells_average,
      int      *cells_median,
      double   *cells_std_deviation,
      int      *cells_min,
      int      *cells_max,
      int      *bound_part_faces_average,
      int      *bound_part_faces_median,
      double   *bound_part_faces_std_deviation,
      int      *bound_part_faces_min,
      int      *bound_part_faces_max,
      int      *bound_part_faces_sum
);

void
PROCF (pdm_part_stat_get, PDM_PART_STAT_GET)
(
const int      *ppartId,
      int      *cells_average,
      int      *cells_median,
      double   *cells_std_deviation,
      int      *cells_min,
      int      *cells_max,
      int      *bound_part_faces_average,
      int      *bound_part_faces_median,
      double   *bound_part_faces_std_deviation,
      int      *bound_part_faces_min,
      int      *bound_part_faces_max,
      int      *bound_part_faces_sum
);


/**
 *
 * \brief Free some array in ppart
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

void
PDM_part_partial_free
(
const  int                ppartId
);

/**
 *
 * \brief Get current meshPat struct
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

PDM_part_t*
PDM_get_meshpart_from_id
(
 int  ppartId
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_part_H__ */
