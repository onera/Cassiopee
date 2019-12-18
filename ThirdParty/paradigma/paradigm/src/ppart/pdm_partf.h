#include "pdmf.h"

  integer, parameter :: PDM_part_SPLIT_PARMETIS = 1
  integer, parameter :: PDM_part_SPLIT_PTSCOTCH = 2
  integer, parameter :: PDM_part_SPLIT_HILBERT = 3

  integer, parameter :: PDM_PART_RENUM_FACE_RANDOM        = 1
  integer, parameter :: PDM_PART_RENUM_FACE_NONE          = 2
  integer, parameter :: PDM_PART_RENUM_FACE_LEXICOGRAPHIC = 3

  integer, parameter :: PDM_part_RENUM_CELL_HILBERT = 1
  integer, parameter :: PDM_part_RENUM_CELL_RANDOM = 2
  integer, parameter :: PDM_part_RENUM_CELL_NONE = 3
  integer, parameter :: PDM_part_RENUM_CELL_CUTHILL = 4

interface

 !================================================================================
 !
 ! \brief Build a initial partitioning
 !
 !  Build a initial partitioning from :
 !      - Cell block distribution with implicit global numbering
 !         (the first cell is the first cell of the first process and
 !          the latest cell is the latest cell of the latest process)
 !      - Face block distribution with implicit global numbering
 !      - Vertex block distribution with implicit global numbering
 !  To repart an existing partition use \ref PDM_part_repart function
 !
 ! \param [out]  ppartId        ppart identifier
 ! \param [in]   pt_comm        MPI Comminicator
 ! \param [in]   method         Choice between (1 for ParMETIS or 2 for PT-Scotch)
 ! \param [in]   nPart          Number of partition to build on this process
 ! \param [in]   dNCell         Number of distributed cells
 ! \param [in]   dNFace         Number of distributed faces
 ! \param [in]   dNVtx          Number of distributed vertices
 ! \param [in]   nFaceGroup     Number of face groups
 ! \param [in]   dCellFaceIdx   Distributed cell face connectivity index or NULL
 !                              (size : dNCell + 1, numbering : 0 to n-1)
 ! \param [in]   dCellFace      Distributed cell face connectivity or NULL
 !                              (size : dFaceVtxIdx[dNCell], numbering : 1 to n)
 ! \param [in]   dCellTag       Cell tag (size : nCell) or NULL
 ! \param [in]   dCellWeight    Cell weight (size : nCell) or NULL
 ! \param [in]   dCellPart      Distributed cell partitioning
 !                              (size = dNCell) or NULL (No partitioning if != NULL)
 ! \param [in]   dFaceCell      Distributed face cell connectivity or NULL
 !                              (size : 2 * dNFace, numbering : 1 to n)
 ! \param [in]   dFaceVtxIdx    Distributed face to vertex connectivity index
 !                              (size : dNFace + 1, numbering : 0 to n-1)
 ! \param [in]   dFaceVtx       Distributed face to vertex connectivity
 !                              (size : dFaceVtxIdx[dNFace], numbering : 1 to n)
 ! \param [in]   dFaceTag       Distributed face tag (size : dNFace)
 !                              or NULL
 ! \param [in]   dVtxCoord      Distributed vertex coordinates
 !                              (size : 3!dNVtx)
 ! \param [in]   dVtxTag        Distributed vertex tag (size : dNVtx) or NULL
 ! \param [in]   dFaceGroupIdx  Index of distributed faces list of each group
 !                              (size = nFaceGroup + 1) or NULL
 ! \param [in]   dFaceGroup     distributed faces list of each group
 !                              (size = dFaceGroup[dFaceGroupIdx[nFaceGroup]],
 !                              numbering : 1 to n)
 !                              or NULL
 !================================================================================

    subroutine pdm_part_create(ppartId, &
                          pt_comm, &
                          method,  &
                          nPart, &
                          dNCell, &
                          dNFace, &
                          dNVtx,&
                          nFaceGroup, &
                          have_dCellFace,&
                          dCellFaceIdx,&
                          dCellFace,&
                          have_dCellTag,&
                          dCellTag, &
                          have_dCellWeight,&
                          dCellWeight,&
                          have_dCellPart,&
                          dCellPart,&
                          have_dFaceCell,&
                          dFaceCell,&
                          dFaceVtxIdx,&
                          dFaceVtx,&
                          have_dFaceTag,&
                          dFaceTag,&
                          dVtxCoord,&
                          have_dVtxTag,&
                          dVtxTag,&
                          dFaceGroupIdx,&
                          dFaceGroup)

    implicit none

    integer                     ::  ppartId
    integer                     ::  pt_comm
    integer                     ::  method
    integer                     ::  nPart
    integer                     ::  dNCell
    integer                     ::  dNFace
    integer                     ::  dNVtx
    integer                     ::  nFaceGroup
    integer                     ::  have_dCellFace
    integer                     ::  dCellFaceIdx(*)
    integer (kind = PDM_g_num_s) :: dCellFace(*)
    integer                     ::  have_dCellTag
    integer                     ::  dCellTag(*)
    integer                     ::  have_dCellWeight
    integer                     ::  dCellWeight(*)
    integer                     ::  have_dCellPart
    integer                     ::  dCellPart(*)
    integer                     ::  have_dFaceCell
    integer (kind = PDM_g_num_s) :: dFaceCell(*)
    integer                     ::  dFaceVtxIdx(*)
    integer (kind = PDM_g_num_s) :: dFaceVtx(*)
    integer                     ::  have_dFaceTag
    integer                     ::  dFaceTag(*)
    double precision            ::  dVtxCoord(*)
    integer                     ::  have_dVtxTag
    integer                     ::  dVtxTag(*)
    integer                     ::  dFaceGroupIdx(*)
    integer (kind = PDM_g_num_s) :: dFaceGroup(*)

  end subroutine


 !==============================================================================
 !
 ! \brief Return a mesh partition dimensions
 !
 ! \param [in]   ppartId            ppart identifier
 ! \param [in]   ipart              Current partition
 ! \param [out]  nCell              Number of cells
 ! \param [out]  nFace              Number of faces
 ! \param [out]  nFacePartBound     Number of partitioning boundary faces
 ! \param [out]  nVtx               Number of vertices
 ! \param [out]  nProc              Number of processus
 ! \param [out]  nTPart             Total number of partitions
 ! \param [out]  sCellFace          Size of cell-face connectivity
 ! \param [out]  sFaceVtx           Size of face-vertex connectivity
 ! \param [out]  sFacePartBound     Size of facePartBound array
 ! \param [out]  sFaceGroup         Size of faceGroup array
 !
 !==============================================================================


   subroutine pdm_part_part_dim_get (ppartId, &
                                  ipart, &
                                  nCell, &
                                  nFace, &
                                  nFacePartBound, &
                                  nVtx, &
                                  nProc, &
                                  nTPart, &
                                  sCellFace, &
                                  sFaceVtx, &
                                  sFaceGroup)
     implicit none

     integer :: ppartId
     integer :: ipart
     integer :: nCell
     integer :: nFace
     integer :: nFacePartBound
     integer :: nVtx
     integer :: nProc
     integer :: nTPart
     integer :: sCellFace
     integer :: sFaceVtx
     integer :: sFaceGroup

   end subroutine pdm_part_part_dim_get


 !==============================================================================
 !
 ! \brief Return a mesh partition
 !
 ! \param [in]   ppartId            ppart identifier
 ! \param [in]   ipart              Current partition
 ! \param [out]  cellTag            Cell tag (size = nCell)
 ! \param [out]  cellFaceIdx        Cell to face connectivity index
 !                                  (size = nCell + 1, numbering : 0 to n-1)
 ! \param [out]  cellFace           Cell to face connectivity
 !                                 (size = cellFaceIdx[nCell] = lCellFace
 !                                  numbering : 1 to n)
 ! \param [out]  cellLNToGN         Cell local numbering to global numbering
 !                                  (size = nCell, numbering : 1 to n)
 ! \param [out]  faceTag            Face tag (size = nFace)
 ! \param [out]  faceCell           Face to cell connectivity
 !                                  (size = 2 * nFace, numbering : 1 to n)
 ! \param [out]  faceVtxIdx         Face to Vertex connectivity index
 !                                  (size = nFace + 1, numbering : 0 to n-1)
 ! \param [out]  faceVtx            Face to Vertex connectivity
 !                                  (size = faceVertexIdx[nFace], numbering : 1 to n)
 ! \param [out]  faceLNToGN         Face local numbering to global numbering
 !                                  (size = nFace, numbering : 1 to n)
 ! \param [out]  facePartBoundProcIdx Partitioning boundary faces block
 !                                    distribution from processus
 !                                    (size = nProc + 1)
 ! \param [out]  facePartBoundPartIdx Partitioning boundary faces block
 !                                    distribution from partition
 !                                    (size = nTPart + 1)
 ! \param [out]  facePartBound      Partitioning boundary faces
 !                                  (size = 4 * nFacePartBound)
 !                                    sorted by processus,
 !                                    sorted by partition in each processus, and
 !                                    sorted by absolute face number in each partition
 !                                  For each face :
 !                                     - Face local number
 !                                       (numbering : 1 to n)
 !                                     - Connected process
 !                                       (numbering : 0 to n-1)
 !                                     - Connected Partition
 !                                       on the connected process
 !                                       (numbering :1 to n)
 !                                     - Connected face local number
 !                                       in the connected partition
 !                                       (numbering :1 to n)
 ! \param [out]  vtxTag             Vertex tag (size = nVertex)
 ! \param [out]  vtx                Vertex coordinates (size = 3 * nVertex)
 ! \param [out]  vtxLNToGN          Vertex local numbering to global numbering
 !                                  (size = nVtx, numbering : 1 to n)
 ! \param [out]  faceGroupIdx       Face group index
 !                                  (size = nFaceGroup + 1, numbering : 1 to n-1)
 ! \param [out]  faceGroup          faces for each group
 !                                  (size = faceGroupIdx[nFaceGroup] = lFaceGroup,
 !                                   numbering : 1 to n)
 ! \param [out]  faceGroupLNToGN    Faces global numbering for each group
 !                                  (size = faceGroupIdx[nFaceGroup] = lFaceGroup,
 !                                  numbering : 1 to n)
 !==============================================================================

 subroutine pdm_part_part_val_get (ppartId, &
                                ipart, &
                                cellTag, &
                                cellFaceIdx, &
                                cellFace, &
                                cellLNToGN, &
                                faceTag, &
                                faceCell, &
                                faceVtxIdx, &
                                faceVtx, &
                                faceLNToGN, &
                                facePartBoundProcIdx, &
                                facePartBoundPartIdx, &
                                facePartBound, &
                                vtxTag, &
                                vtx, &
                                vtxLNToGN, &
                                faceGroupIdx, &
                                faceGroup, &
                                faceGroupLNToGN)

   implicit none

   integer                       :: ppartId
   integer                       :: ipart
   integer                       :: cellTag(*)
   integer                       :: cellFaceIdx(*)
   integer                       :: cellFace(*)
   integer (kind = PDM_g_num_s)  :: cellLNToGN(*)
   integer                       :: faceTag(*)
   integer                       :: faceCell(*)
   integer                       :: faceVtxIdx(*)
   integer                       :: faceVtx(*)
   integer (kind = PDM_g_num_s)  :: faceLNToGN(*)
   integer                       :: facePartBoundProcIdx(*)
   integer                       :: facePartBoundPartIdx(*)
   integer                       :: facePartBound(*)
   integer                       :: vtxTag(*)
   double precision              :: vtx(*)
   integer (kind = PDM_g_num_s)  :: vtxLNToGN(*)
   integer                       :: faceGroupIdx(*)
   integer                       :: faceGroup(*)
   integer (kind = PDM_g_num_s)  :: faceGroupLNToGN(*)

 end subroutine pdm_part_part_val_get

 !==============================================================================
 !
 ! \brief Free ppart
 !
 ! \param [in]   ppartId        ppart identifier
 !
 !==============================================================================

 subroutine pdm_part_free (ppartId)

   implicit none

   integer                       :: ppartId

 end subroutine pdm_part_free

 !==============================================================================
 !
 ! \brief Return times
 !
 ! \param [in]   ppartId     ppart identifier
 ! \param [out]  elapsed     elapsed times (size = 4)
 ! \param [out]  cpu         cpu times (size = 4)
 ! \param [out]  cpu_user    user cpu times (size = 4)
 ! \param [out]  cpu_sys     system cpu times (size = 4)
 !
 !==============================================================================

 subroutine pdm_part_time_get (ppartId, &
                            elapsed, &
                            cpu, &
                            cpu_user, &
                            cpu_sys)
   implicit none

   integer           :: ppartId
   double precision  :: elapsed
   double precision  :: cpu
   double precision  :: cpu_user
   double precision  :: cpu_sys

 end subroutine pdm_part_time_get

 !==============================================================================
 !
 ! \brief Return statistic
 !
 ! \param [in]   ppartId                        ppart identifier
 ! \param [out]  cells_average                  average of cells number
 ! \param [out]  cells_median                   median of cells number
 ! \param [out]  cells_std_deviation            standard deviation of cells number
 ! \param [out]  cells_min                      minimum of cells nummber
 ! \param [out]  cells_max                      maximum of cells nummber
 ! \param [out]  bound_part_faces_average       average of partitioning boundary faces
 ! \param [out]  bound_part_faces_median        median of partitioning boundary faces
 ! \param [out]  bound_part_faces_std_deviation standard deviation of partitioning boundary faces
 ! \param [out]  bound_part_faces_min           minimum of partitioning boundary faces
 ! \param [out]  bound_part_faces_max           maximum of partitioning boundary faces
 !
 !==============================================================================

 subroutine pdm_part_stat_get (ppartId,  &
                            cells_average, &
                            cells_median, &
                            cells_std_deviation, &
                            cells_min,  &
                            cells_max, &
                            bound_part_faces_average, &
                            bound_part_faces_median,  &
                            bound_part_faces_std_deviation, &
                            bound_part_faces_min, &
                            bound_part_faces_max, &
                            bound_part_faces_sum)

   implicit none

   integer      :: ppartId
   integer      :: cells_average
   integer      :: cells_median
   double precision   :: cells_std_deviation
   integer      :: cells_min
   integer      :: cells_max
   integer      :: bound_part_faces_average
   integer      :: bound_part_faces_median
   double precision   :: bound_part_faces_std_deviation
   integer      :: bound_part_faces_min
   integer      :: bound_part_faces_max
   integer      :: bound_part_faces_sum

 end subroutine pdm_part_stat_get


end interface
