module mod_pdm_part

  use mod_pdm

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
  
  interface pdm_part_create     ; module procedure  &
    pdm_part_create_
  end interface

  interface pdm_part_coarse_mesh_create     ; module procedure  &
    pdm_part_coarse_mesh_create_
  end interface

  interface pdm_part_renum_method_cell_idx_get      ; module procedure  &
    pdm_part_renum_method_cell_idx_get_
  end interface 

  interface pdm_part_renum_method_face_idx_get      ; module procedure  &
    pdm_part_renum_method_face_idx_get_
  end interface  

  interface pdm_part_renum_method_cell_name_get      ; module procedure  &
    pdm_part_renum_method_cell_name_get_
  end interface  

  interface pdm_part_renum_method_face_name_get      ; module procedure  &
    pdm_part_renum_method_face_name_get_
  end interface  
  
interface

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
                                  nProc, &
                                  nTPart, &
                                  nVtx, &
                                  sCellFace, &
                                  sFaceVtx, &
                                  sFaceGroup, &
                                  nFaceGroup)
     use mod_pdm

     implicit none
     
     integer :: ppartId
     integer :: ipart
     integer :: nCell
     integer :: nFace
     integer :: nFacePartBound
     integer :: nProc
     integer :: nTPart
     integer :: nVtx
     integer :: sCellFace
     integer :: sFaceVtx
     integer :: sFaceGroup
     integer :: nFaceGroup

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

   use mod_pdm

   implicit none

   integer                       :: ppartId
   integer                       :: ipart
   integer                       :: cellTag(*)
   integer                       :: cellFaceIdx(*)
   integer                       :: cellFace(*)
   integer (kind = PDM_g_num_s) :: cellLNToGN(*)
   integer                       :: faceTag(*)
   integer                       :: faceCell(*)
   integer                       :: faceVtxIdx(*)
   integer                       :: faceVtx(*)
   integer (kind = PDM_g_num_s) :: faceLNToGN(*)
   integer                       :: facePartBoundProcIdx(*)
   integer                       :: facePartBoundPartIdx(*)
   integer                       :: facePartBound(*)
   integer                       :: vtxTag(*)
   double precision              :: vtx(*)
   integer (kind = PDM_g_num_s) :: vtxLNToGN(*)
   integer                       :: faceGroupIdx(*)
   integer                       :: faceGroup(*)
   integer (kind = PDM_g_num_s) :: faceGroupLNToGN(*)
 
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
   double precision  :: elapsed(4)
   double precision  :: cpu(4)
   double precision  :: cpu_user(4)
   double precision  :: cpu_sys(4)

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
      
   use mod_pdm
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

 
 !================================================================================
 !
 ! \brief Get the number of renumbering cell methods 
 ! 
 ! \param [out]    Number of methods
 !
 !================================================================================

 subroutine pdm_part_n_renum_method_cell_get (n_method)
   use mod_pdm
   implicit none
   integer      :: n_method

 end subroutine pdm_part_n_renum_method_cell_get

 
 
 !================================================================================
 !
 ! \brief Get the number of renumbering cell methods 
 ! 
 ! \param [out]    Number of methods
 !
 !================================================================================

 subroutine pdm_part_n_renum_method_face_get (n_method)
   use mod_pdm
   implicit none
   integer      :: n_method

 end subroutine pdm_part_n_renum_method_face_get
 
end interface

private :: pdm_part_create_ ,&
           pdm_part_coarse_mesh_create_ ,&
           pdm_part_renum_method_cell_idx_get_ ,&
           pdm_part_renum_method_face_idx_get_ ,&
           pdm_part_renum_method_cell_name_get_ ,&
           pdm_part_renum_method_face_name_get_

contains 

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

   subroutine pdm_part_create_(ppartId, &
                          pt_comm, &
                          split_method,  &
                          renum_cell_method, &
                          renum_face_method, &
                          nPropertyCell, &
                          renum_properties_cell, &
                          nPropertyFace, &
                          renum_properties_face, &
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

    use mod_pdm

    implicit none

    integer                     ::  ppartId
    integer                     ::  pt_comm
    integer                     ::  split_method
    character (len=*)           ::  renum_cell_method
    character (len=*)           ::  renum_face_method
    integer                     ::  nPropertyCell
    integer                     ::  renum_properties_cell(*)
    integer                     ::  nPropertyFace
    integer                     ::  renum_properties_face(*)
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
    integer                     :: l_renum_cell_method
    integer                     :: l_renum_face_method
    
    l_renum_cell_method = len(renum_cell_method)
    l_renum_face_method = len(renum_face_method)

    call pdm_part_create_cf (ppartId, &
                             pt_comm, &
                             split_method, &
                             renum_cell_method, &
                             l_renum_cell_method, & 
                             renum_face_method, &
                             l_renum_face_method, & 
                             nPropertyCell, &
                             renum_properties_cell, &
                             nPropertyFace, &
                             renum_properties_face, &
                             nPart, &
                             dNCell, &
                             dNFace, &
                             dNVtx, &
                             nFaceGroup, &
                             have_dCellFace, &
                             dCellFaceIdx, &
                             dCellFace, &
                             have_dCellTag, &
                             dCellTag, &
                             have_dCellWeight, &
                             dCellWeight, &
                             have_dCellPart, &
                             dCellPart, &
                             have_dFaceCell, &
                             dFaceCell, &
                             dFaceVtxIdx, &
                             dFaceVtx, &
                             have_dFaceTag, &
                             dFaceTag, &
                             dVtxCoord, &
                             have_dVtxTag, &
                             dVtxTag, &
                             dFaceGroupIdx, &
                             dFaceGroup)

  end subroutine pdm_part_create_


 !================================================================================
 !
 ! \brief Return an initialized coarse mesh object
 !
 ! \param [out]  cmId              Coarse mesh identifier
 ! 
 ! \param [in]   pt_comm           Communicator
 ! \param [in]   method            Choice between (1 for ParMETIS or 2 for PT-Scotch)
 ! \param [in]   nPart             Number of partitions
 ! \param [in]   nTPart            Total number of partitions
 ! \param [in]   have_cellTag      Presence d'un tableau de tags pour les cellules
 ! \param [in]   have_faceTag      Presence d'un tableau de tags pour les faces
 ! \param [in]   have_vtxTag       Presence d'un tableau de tags pour les sommets
 ! \param [in]   have_cellWeight   Presence d'un tableau de poids pour les cellules
 ! \param [in]   have_faceWeight   Presence d'un tableau de poids pour les faces
 ! \param [in]   have_faceGroup    Presence des tableaux de groupes de faces
 !================================================================================

    
subroutine pdm_part_coarse_mesh_create_ (cmId, &
                                        comm, &        
                                        method, &
                                        nPart, &
                                        nTPart, &
                                        nFaceGroup,&
                                        have_cellTag,&
                                        have_faceTag,&
                                        have_vtxTag,&
                                        have_cellWeight, &
                                        have_faceWeight, &
                                        have_faceGroup)

    use mod_pdm

    implicit none

    integer                     ::  cmId
    integer                     ::  comm
    character (len=*)           ::  method
    integer                     ::  nPart
    integer                     ::  nTPart
    integer                     ::  nFaceGroup
    integer                     ::  have_cellTag
    integer                     ::  have_faceTag
    integer                     ::  have_vtxTag
    integer                     ::  have_cellWeight
    integer                     ::  have_faceWeight
    integer                     ::  have_faceGroup

    integer                     :: l_method
    
    l_method = len(method)
  

    call pdm_part_coarse_mesh_create_cf (cmId, &
                                         comm, &        
                                         method, &
                                         l_method, &
                                         nPart, &
                                         nTPart, &
                                         nFaceGroup,&
                                         have_cellTag,&
                                         have_faceTag,&
                                         have_vtxTag,&
                                         have_cellWeight, &
                                         have_faceWeight, &
                                         have_faceGroup)

    
 end subroutine pdm_part_coarse_mesh_create_

  
 !================================================================================
 !
 ! \brief Get index of a renumbering face method
 ! 
 ! \param [in]       name   Name of the method
 ! \param [in, out]  idx    Index of method -1 otherwise
 ! 
 !================================================================================

  subroutine pdm_part_renum_method_face_idx_get_ (name, &
                                                  idx)

     use mod_pdm

    implicit none

    character (len=*) :: name
    integer           :: idx
    
    integer           :: l_name
    
    l_name = len(name)

    call pdm_part_renum_method_face_idx_get_cf (name, l_name, idx)

  end subroutine pdm_part_renum_method_face_idx_get_


 !================================================================================
 !
 ! \brief Get index of a renumbering cell method
 ! 
 ! \param [in]       name   Name of the method
 ! \param [in, out]  idx    Index of method -1 otherwise
 ! 
 !================================================================================

 subroutine pdm_part_renum_method_cell_idx_get_ (name, &
                                                 idx)

   use mod_pdm

   implicit none

   character (len=*) :: name
   integer           :: idx
   
   integer           :: l_name

   l_name = len(name)
   
   call pdm_part_renum_method_cell_idx_get_cf (name, l_name, idx)

 end subroutine pdm_part_renum_method_cell_idx_get_

 !================================================================================
 !
 ! \brief Get name of the face renumbering method 
 ! 
 ! \param [in]  idx     Index of the method
 ! \param [in, out]     Name  of the method, '' otherwize
 ! 
 !================================================================================


subroutine pdm_part_renum_method_face_name_get_ (idx, &
                                                name)

   use mod_pdm
   implicit none

   character (len = *) :: name
   integer  :: idx

   integer  :: l_name
   l_name = len(name)

   call pdm_part_renum_method_face_name_get_cf (name, l_name, idx)


end subroutine pdm_part_renum_method_face_name_get_
  

 !================================================================================
 !
 ! \brief Get name of the face renumbering method 
 ! 
 ! \param [in]  idx     Index of the method
 ! \param [in, out]     Name  of the method, '' otherwize
 ! 
 !================================================================================


subroutine pdm_part_renum_method_cell_name_get_ (idx, &
                                                name)

   use mod_pdm
   implicit none

   character (len = *) :: name
   integer  :: idx

   integer  :: l_name
   l_name = len(name)

   call pdm_part_renum_method_cell_name_get_cf (name, l_name, idx)


end subroutine pdm_part_renum_method_cell_name_get_

end module mod_pdm_part
