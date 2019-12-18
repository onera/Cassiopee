#ifndef __PDM_PART_PRIV_H__
#define __PDM_PART_PRIV_H__

#include <stdlib.h>

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_timer.h"
#include "pdm_part.h"
#include "pdm_mpi.h"

/*============================================================================
 * Type definitions
 *============================================================================*/


/**
 * \struct _subpartlayout_t
 * \brief  Partition object
 *
 * _subpartlayout_t define a mesh partition layouts base on sudomaine
 *
 */

typedef struct  _subpartlayout_t {

  int           nSdom;                  /*!< Number of subDomain                     */
  int           nFaceInt;               /*!< Number of Interior face                 */
  int           nFaceExt;               /*!< Number of Exterior face                 */

  /* Idx array of displacement */
  int*          cellTileIdx;           /*!< Cell Tile Index     (Size = nSdom + 1)   */
  int*          faceTileIdx;           /*!< Face Tile Index     (Size = nSdom + 1)   */
  int*          faceBndTileIdx;        /*!< Face Bnd Tile Index (Size = nSdom + 1)   */

  /* Idx array of displacement */
  int*          maskTileIdx;           /*!< Mask Tile Index   (Size = nSdom + 1)     */
  int*          cellVectTileIdx;       /*!< Cell Tile Index   (Size = nSdom + 1)     */
  int*          maskTileN;             /*!< Mask Tile number  (Size = nSdom + 1)     */
  int*          cellVectTileN;         /*!< Cell Tile number  (Size = nSdom + 1)     */
  int*          maskTile;              /*!< Mask Tile number                         */


} _subpartlayout_t;


/**
 * \struct _part_t
 * \brief  Partition object
 *
 * _part_t define a mesh partition
 *
 */

typedef struct  _part_t {
  int           nVtx;               /*!< Number of vertices                   */
  int           nCell;              /*!< Number of cells                      */
  int           nFace;              /*!< Number of faces                      */
  int           nFacePartBound;     /*!< Number of partitioning boundary faces*/
  int           nFaceGroup;         /*!< Number of boundary faces             */

  int          *cellFaceIdx;        /*!< Cell face connectivity index
                                      (size = nCell + 1)                      */
  PDM_g_num_t *gCellFace;          /*!< Global numbering cell face connectivity
                                      (size = cellFaceIdx[nCell])             */
  int          *cellFace;           /*!< Cell face connectivity
                                      (size = cellFaceIdx[nCell])             */
  PDM_g_num_t *cellLNToGN;         /*!< Local to global cell numbering
                                      (size = nCell)                          */
  int          *cellTag;            /*!< Cell tag
                                      (size = nCell)                          */

  int          *faceCell;           /*!< Face cell connectivity
                                      (size = 2 * nFace)                      */
  int          *faceVtxIdx;         /*!< Face vertex connectivity index
                                      (size = nFace + 1)                      */
  PDM_g_num_t *gFaceVtx;           /*!< Global numbering face vtx connectivity
                                      (size = faceVtxIdx[nFace])              */
  int          *faceVtx;            /*!< Face vertex connectivity
                                      (size = faceVtxIdx[nFace])              */
  PDM_g_num_t *faceLNToGN ;         /*!< Local to global cell numbering
                                      (size = nFace)                          */
  int          *faceTag;            /*!< Tag face connectivity index
                                      (size = nFace)                          */

  int          *facePartBoundProcIdx;   /*!< Partitioning boundary bloc distribution
                                      (size = nRanks + 1)                     */
  int          *facePartBoundPartIdx;   /*!< Partitioning boundary bloc distribution
                                      (size = nPartTotal + 1)                 */
  int          *facePartBound;      /*!< Partitioning boundary faces sorted by
                                         proc, sorted by part in the proc, and
                                         sorted by absolute face number in the part
                                         For each face :
                                          - Face local number
                                          - Connected process
                                          - Connected Partition
                                            on the connected process
                                          - Connected face local number
                                            in the connected partition
                                      (size = 4* nFacePartBound)              */

  int          *faceGroupIdx;       /*!< Face group index
                                      (size = nFaceGroup + 1)                     */
  int          *faceGroup;          /*!< Face group index
                                      (size = faceGroupIdx[nFaceGroup])       */

  PDM_g_num_t *faceGroupLNToGN;    /*!< Local to global boundary face numbering
                                      (size = faceGroupIdx[nFaceGroup])           */

  double       *vtx;                /*!< Vertex coordinates
                                      (size = 3 * nVtx)                           */
  PDM_g_num_t *vtxLNToGN;          /*!< Local to global vertex numbering
                                      (size = nVtx)                           */
  int          *vtxTag;             /*!< Tag vertex
                                      (size = nVtx)                           */


  const int          *cellWeight;             /*!< Cell weight - For coarse mesh
                                            (size = nCel)                           */
  const int          *faceWeight;             /*!< Face weight - For coarse mesh
                                            (size = nFac)                           */

  int          *cellColor;             /*!< Cell color - For cache blocking
                                            (size = nCel)                           */
  int          *faceColor;             /*!< Face color - For cache blocking
                                            (size = nFac)                           */
  int          *threadColor;             /*!< Thread color - For cache blocking
                                            (size = nThread)                        */
  int          *hyperPlaneColor;         /*!< Thread color - For cache blocking
                                            (size = nThread)                         */

  int          *newToOldOrderCell;   /*!< Cell reordering
                                         (size = nCel)                           */
  int          *newToOldOrderFace;   /*!< Face reordering
                                            (size = nFac)                        */

  _subpartlayout_t *subpartlayout;    /*!< Layouts of subdomain                     */

} _part_t;

/**
 * \struct _PDM_part_t
 * \brief   PPART object
 *
 * _part_t define a parallel mesh partitioning
 *
 */

typedef struct _PDM_part_t {

  /* Local dimensions */

  int                 dNVtx;         /*!< Number of distributed vertices      */
  int                 dNCell;        /*!< Number of distributed cells         */
  int                 dNFace;        /*!< Number of distributed faces         */
  int                 nFaceGroup;    /*!< Number of boundaries                */

  /* Cell definitions */

  const int          *_dCellFaceIdx;  /*!< Cell-face connectivity of distributed
                                       cells (size = dcellFaceIdx[dNCell],
                                       shared array) (computed)               */
  const PDM_g_num_t *_dCellFace;    /*!< Tag of distributed cells
                                       (size = dNCell, shared array)          */
  const int          *_dCellTag;     /*!< Tag of distributed cells
                                       (size = dNCell, shared array)          */
  const int          *_dCellWeight;  /*!< Weight of distributed cells
                                       (size = dNCell, shared array)          */
  const int          *_dCellPart;    /*!< Partitioning of distributed cells
                                       (size = dNCell, shared array)          */
  int                *dCellFaceIdx;  /*!< Cell-face connectivity index of
                                       distributed cells
                                       (size = dNCell + 1)
                                        computed                             */
  PDM_g_num_t       *dCellFace;     /*!< Cell-face connectivity of distributed
                                       cells (size = dcellFaceIdx[dNCell],
                                       computed                               */
  PDM_g_num_t       *dCellProc;     /*!< Initial cell distribution on processes
                                       (size = nRank + 1) (computed)          */

  /* Face definitions */

  const int          *_dFaceTag;     /*!< Tag of distributed face
                                       (size = dNFace, shared array)          */
  const PDM_g_num_t *_dFaceCell;    /*!< Face-cell connectivity of distributed
                                       faces (size = 2 * dNFace, shared array)
                                       if iface is a boundary face,
                                       _dFaceCell[2*iface + 1] = 0           */
  const int          *_dFaceVtxIdx;  /*!< Face-vertex connectivity index of
                                       distributed faces (size = dNFace + 1,
                                       shared array)                          */
  const PDM_g_num_t *_dFaceVtx;     /*!< Face-vertex connectivity of
                                       distributed faces
                                       (size = dFaceVtxIdx[dNFace],shared array)
                                     */
  PDM_g_num_t       *dFaceProc;     /*!< Initial face distribution on processes
                                       (size = nRank + 1)                     */

  PDM_g_num_t       *dFaceCell;    /*!< Face-cell connectivity of distributed
                                       faces (size = 2 * dNFace,) computed
                                       if iface is a boundary face,
                                       _dFaceCell[2*iface + 1] = 0           */
  /* Vertex definitions */

  const double       *_dVtxCoord;    /*!< Coordinates of ditributed vertices
                                       (size = 3 * dNVtx, shared array)       */
  const int          *_dVtxTag;      /*!< Tag of distributed vertices
                                       (size = dNVtx, shared array)           */
  PDM_g_num_t       *dVtxProc;      /*!< Initial vertex distribution
                                       (size = nrank + 1                      */

  /* Face group */

  const int          *_dFaceGroupIdx; /*!< Index of distributed faces list of
                                        each boundary (size = nBound + 1)
                                        or NULL                               */
  const PDM_g_num_t *_dFaceGroup    ;/*!< Distributed faces list of each
                                       boundary (size = dfaceBoundIdx[nBound])
                                       or NULL                                */

  /* Partitioning boundary faces */

  int          *dPartBound;          /*!< Partitioning boundary faces definition
                                       For each face :
                                       (rank1, part1, lFace1, rank2, part2, lFace2)
                                       or -1 size = 6 * dNface  */

  /* Dual graph */

  PDM_g_num_t *dDualGraphIdx;      /*!< Dual graph index
                                      (size = dNCell + 1)                     */
  PDM_g_num_t *dDualGraph;         /*!< Dual graph
                                      (size = dualGraphIdx[dNCell])           */

  PDM_timer_t *timer;             /*!< Timer */

  double times_elapsed [4];          /*!< Elapsed times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */

  double times_cpu[4];             /*!< CPU times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */

  double times_cpu_u[4];           /*!< User CPU times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */

  double times_cpu_s[4];          /*!< Systeme CPU times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */

  /* Communicator */

  PDM_MPI_Comm      comm;               /*!< Communicator */

  /* Partitions */

  PDM_part_split_t split_method;             /*!< Partitioning method */

  int renum_face_method;                     /*!< Renumbering face method */

  int renum_cell_method;                     /*!< Renumbering cell method */

  int  nPropertyCell;                         /*!< Size of cells properties      */
  int  nPropertyFace;                         /*!< Size of faces properties      */
  const int* renum_properties_cell;           /*!< Renumbering cells properties  */
  const int* renum_properties_face;           /*!< Renumbering faces properties  */


  int          nPart;               /*!< Number of partitions to define
                                      on this process */
  int          tNPart;              /*!< Total number of partitions           */
  int          mNPart;              /*!< Maximum number of partitions to define
                                      on one process */
  int         *dPartProc;           /*!< Initial cell distribution on processes
                                      (size = nRank + 1)                      */
  int         *gPartTolProcPart;    /*!< For each lobal part number :
                                          - process storing this partition
                                          - local number partition on this
                                            partition
                                         (size = 2*tNpart)                    */
 _part_t     **meshParts;           /*!< Partitions built ont this process
                                      (size = npart)                          */

} _PDM_part_t;

/**
 *
 * \brief Return an initialized part object
 *
 */

static inline _part_t *
_part_create
(
void
 )
{
  _part_t *part = (_part_t *) malloc(sizeof(_part_t));

  part->nVtx = 0;
  part->nCell = 0;
  part->nFace = 0;
  part->nFacePartBound = 0;
  part->cellFaceIdx = NULL;
  part->gCellFace = NULL;
  part->cellFace = NULL;
  part->cellLNToGN = NULL;
  part->cellTag = NULL;
  part->faceCell = NULL;
  part->faceVtxIdx = NULL;
  part->gFaceVtx = NULL;
  part->faceVtx = NULL;
  part->faceLNToGN = NULL;
  part->faceTag = NULL;
  part->facePartBoundProcIdx = NULL;
  part->facePartBoundPartIdx = NULL;
  part->facePartBound = NULL;
  part->faceGroupIdx = NULL;
  part->faceGroup = NULL;
  part->faceGroupLNToGN = NULL;
  part->vtx = NULL;
  part->vtxLNToGN = NULL;
  part->vtxTag = NULL;
  part->cellColor = NULL;
  part->faceColor = NULL;
  part->threadColor = NULL;
  part->hyperPlaneColor = NULL;
  part->newToOldOrderCell = NULL;
  part->newToOldOrderFace = NULL;
  part->subpartlayout = NULL;
  return part;
}

#endif
