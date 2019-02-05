
cdef extern from "pdm_dmesh_nodal.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure 
    ctypedef struct PDM_DMesh_nodal_t:
      pass
      
    ctypedef enum PDM_Mesh_nodal_elt_t: 
      PDM_MESH_NODAL_POINT    = 1   
      PDM_MESH_NODAL_BAR2     = 2  
      PDM_MESH_NODAL_TRIA3    = 3   
      PDM_MESH_NODAL_QUAD4    = 4   
      PDM_MESH_NODAL_POLY_2D  = 5     
      PDM_MESH_NODAL_TETRA4   = 6    
      PDM_MESH_NODAL_PYRAMID5 = 7      
      PDM_MESH_NODAL_PRISM6   = 8    
      PDM_MESH_NODAL_HEXA8    = 9   
      PDM_MESH_NODAL_POLY_3D  = 10 
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function 
    int PDM_DMesh_nodal_create(PDM_MPI_Comm comm, PDM_g_num_t nCel, PDM_g_num_t nVtx)
    void PDM_DMesh_nodal_free(int handle)
    
    void PDM_DMesh_nodal_coord_set(int handle, int n_vtx, double* coords)

    PDM_g_num_t *PDM_DMesh_nodal_distrib_vtx_get(int handle)
    PDM_g_num_t *PDM_DMesh_nodal_distrib_section_get(int handle, int id_section)
    
    int                  PDM_DMesh_nodal_n_vtx_get(int handle)
    int                  PDM_DMesh_nodal_n_sections_get(int handle)
    int*                 PDM_DMesh_nodal_sections_id_get(int handle)
    PDM_Mesh_nodal_elt_t PDM_DMesh_nodal_section_type_get(int handle, int id_section)
    double*              PDM_DMesh_nodal_vtx_get(int handle)
    
    int                  PDM_DMesh_nodal_section_add(int handle, PDM_Mesh_nodal_elt_t t_elt)
    void                 PDM_DMesh_nodal_section_std_set(int          handle, 
                                                        int          id_section,
                                                        int          n_elmts, 
                                                        PDM_g_num_t* connec)
    
    PDM_g_num_t* PDM_DMesh_nodal_section_std_get(int          handle, int id_section)
    int PDM_DMesh_nodal_section_n_elt_get(int          handle, int id_section)
    
    void PDM_DMesh_nodal_section_poly2d_set(int handle, int id_section, PDM_l_num_t n_elt, 
                                            PDM_l_num_t* connec_idx, 
                                            PDM_g_num_t   *connec)
    
    PDM_g_num_t PDM_DMesh_nodal_total_n_cell_get(int handle)
    PDM_g_num_t PDM_DMesh_nodal_total_n_face_get(int handle)
    PDM_g_num_t PDM_DMesh_nodal_total_n_vtx_get(int handle)
    void PDM_DMesh_nodal_cell_face_compute(int handle)
    int PDM_DMesh_nodal_cell_face_get(int handle, int** cell_face_idx, PDM_g_num_t **cell_face)
    int PDM_DMesh_nodal_face_cell_get(int handle, PDM_g_num_t** face_cell)
    int PDM_DMesh_nodal_face_vtx_get(int handle, int** dface_vtx_idx, PDM_g_num_t **dface_vtx)
    
    PDM_g_num_t* PDM_DMesh_nodal_distrib_cell_get(int handle)
    PDM_g_num_t* PDM_DMesh_nodal_distrib_face_get(int handle)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

cdef extern from "pdm_elt_parent_find.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure     
    void PDM_elt_parent_find_from_distrib(PDM_g_num_t  *elt_distrib,
                                          int          *elt_def_idx,
                                          PDM_g_num_t  *elt_def,
                                          PDM_g_num_t  *elt_to_find_distrib,
                                          int          *elt_to_find_def_idx,
                                          PDM_g_num_t  *elt_to_find_def,
                                          PDM_MPI_Comm  comm,     
                                          PDM_g_num_t  *parent)
    
    void PDM_elt_parent_find(int           dnelt,
                             int          *elt_def_idx,
                             PDM_g_num_t  *elt_def,
                             int           dnelt_to_find,
                             int          *elt_to_find_def_idx,
                             PDM_g_num_t  *elt_to_find_def,
                             PDM_MPI_Comm  comm,     
                             PDM_g_num_t  *parent)
    
cdef extern from "pdm_distrib.h":

    void PDM_distrib_compute(int           dnelt, 
                             PDM_g_num_t  *elt_distrib,   
                             int           offset, 
                             PDM_MPI_Comm  comm)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ------------------------------------------------------------------
cdef class DistributedMeshNodal:
    """
       DistributedMeshNodal: Interface to build face from Element->Vtx connectivity
    """
    # ************************************************************************
    # > Class attributes
    cdef int idmesh
    cdef int Size
    cdef int Rank
    # ************************************************************************
    # ------------------------------------------------------------------------
    def __init__(self, MPI.Comm    comm, 
                       PDM_g_num_t nVtx, 
                       PDM_g_num_t nCel):
        """
        TODOUX
        """
        # ************************************************************************
        # > Declaration
        # cdef int      nElts
        # cdef int      idx
        # # > Numpy array
        # cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='fortran'] partLNToGN
        # ************************************************************************
        
        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        self.Rank = comm.Get_rank()
        self.Size = comm.Get_size()
        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        
        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Convert mpi4py -> PDM_MPI
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi
        cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        
        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        self.idmesh = PDM_DMesh_nodal_create(PDMC, nVtx, nCel)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def SetCoordinnates(self, NPY.ndarray[NPY.double_t  , mode='c', ndim=1] dVtxCoord):
        """
        """
        # ************************************************************************
        # > Declaration
        cdef int nVtx
        # ************************************************************************

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        nVtx = dVtxCoord.shape[0]
        PDM_DMesh_nodal_coord_set(self.idmesh, nVtx, <double *> dVtxCoord.data)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        
    # ------------------------------------------------------------------------
    def SetSections(self, list ElmtList, 
                          NPY.ndarray[NPY.int32_t, mode='c', ndim=1] ElmtsTyp, 
                          NPY.ndarray[NPY.int32_t, mode='c', ndim=1] nElemts):
        """
           TODO : Split function as PDM
        """
        # ************************************************************************
        # > Declaration
        cdef int nVtx
        cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='fortran'] Connect
        # ************************************************************************
        
        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Panic assert
        assert(len(ElmtList) == ElmtsTyp.shape[0])
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        for iElmt, Connect in enumerate(ElmtList):
          id_section = PDM_DMesh_nodal_section_add(self.idmesh, <PDM_Mesh_nodal_elt_t> ElmtsTyp[iElmt])
          PDM_DMesh_nodal_section_std_set(self.idmesh, id_section, nElemts[iElmt], <PDM_g_num_t *> Connect.data)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def Compute(self):
        """
        """
        # ************************************************************************
        # > Declaration
        # ************************************************************************
        PDM_DMesh_nodal_cell_face_compute(self.idmesh)      
        
    # ------------------------------------------------------------------------
    def getFaceCell(self):
        """
        """
        # ************************************************************************
        # > Declaration
        cdef PDM_g_num_t *faceCell
        cdef PDM_g_num_t nFace
        cdef PDM_g_num_t dNface
        cdef NPY.npy_intp dim
        # ************************************************************************
        
        # > Get Size
        nFace  = PDM_DMesh_nodal_total_n_face_get(self.idmesh)
        
        # > Get array        
        dNface = PDM_DMesh_nodal_face_cell_get(self.idmesh, &faceCell)
        
        # > Build numpy capsule
        dim = <NPY.npy_intp> 2 * dNface
        npFaceCell = NPY.PyArray_SimpleNewFromData(1,
                                                   &dim,
                                                   PDM_G_NUM_NPY_INT,
                                                   <void *> faceCell)
        
        return {'sFace'        : nFace, 
                 'dNFace'      : dNface,
                 'npdFaceCell' : npFaceCell}
    
    # ------------------------------------------------------------------------
    def getCellFace(self):
        """
        """
        # ************************************************************************
        # > Declaration
        cdef PDM_g_num_t *CellFace
        cdef PDM_l_num_t *CellFaceIdx
        cdef PDM_g_num_t nFace
        cdef PDM_g_num_t dNface
        cdef NPY.npy_intp dim
        # ************************************************************************
        
        # > Get Size
        # nCell  = PDM_DMesh_nodal_total_n_cell_get(self.idmesh)
        nCell  = 0
        
        # > Get array   
        dNCell = PDM_DMesh_nodal_cell_face_get(self.idmesh,  &CellFaceIdx, &CellFace)
        
        # > Build numpy capsule
        dim = <NPY.npy_intp> dNCell + 1
        npCellFaceIdx = NPY.PyArray_SimpleNewFromData(1,
                                                      &dim,
                                                      NPY.NPY_INT32,
                                                      <void *> CellFaceIdx)
        
        dim = <NPY.npy_intp> npCellFaceIdx[npCellFaceIdx.shape[0]-1]
        npCellFace = NPY.PyArray_SimpleNewFromData(1,
                                                   &dim,
                                                   PDM_G_NUM_NPY_INT,
                                                   <void *> CellFace)
        
        return {'sCell'           : nCell, 
                'dNCell'          : dNCell,
                'npdCellFaceIdx'  : npCellFaceIdx,
                'npdCellFace'     : npCellFace}
                
    # ------------------------------------------------------------------------
    def getFaceVtx(self):
        """
        """
        # ************************************************************************
        # > Declaration
        cdef PDM_g_num_t *faceVtx
        cdef PDM_l_num_t *faceVtxIdx
        cdef PDM_g_num_t nFace
        cdef PDM_g_num_t dNface
        cdef NPY.npy_intp dim
        # ************************************************************************
        
        # > Get Size
        nFace  = PDM_DMesh_nodal_total_n_face_get(self.idmesh)
        
        # > Get array        
        dNface = PDM_DMesh_nodal_face_vtx_get(self.idmesh, &faceVtxIdx, &faceVtx)
        
        # > Build numpy capsule
        dim = <NPY.npy_intp> dNface + 1
        npfaceVtxIdx = NPY.PyArray_SimpleNewFromData(1,
                                                   &dim,
                                                   PDM_G_NUM_NPY_INT,
                                                   <void *> faceVtxIdx)
        
        dim = <NPY.npy_intp> npfaceVtxIdx[npfaceVtxIdx.shape[0]-1]
        npfaceVtx = NPY.PyArray_SimpleNewFromData(1,
                                                   &dim,
                                                   PDM_G_NUM_NPY_INT,
                                                   <void *> faceVtx)
        
        return {'sFace'         : nFace, 
                 'npdFaceVtxIdx' : npfaceVtxIdx,
                 'npdFaceVtx'    : npfaceVtx}
                 
    # ------------------------------------------------------------------------
    def getDistribFace(self):
        """
        """
        # ************************************************************************
        # > Declaration
        cdef PDM_g_num_t *faceDistrib
        cdef NPY.npy_intp dim
        # ************************************************************************
        
        # > Get array        
        faceDistrib = PDM_DMesh_nodal_distrib_face_get(self.idmesh)
        
        # > Build numpy capsule
        dim = <NPY.npy_intp> self.Rank + 1
        npDistribFace = NPY.PyArray_SimpleNewFromData(1,
                                                      &dim,
                                                      PDM_G_NUM_NPY_INT,
                                                      <void *> faceDistrib)
        
        return npDistribFace

    # ------------------------------------------------------------------------
    def __dealloc__(self):
      """
         Use the free method of PDM Lib
      """
      # ************************************************************************
      # > Declaration
      # ************************************************************************

      # > Free Ppart Structure
      print 'PDM_DMesh_nodal_free'
      PDM_DMesh_nodal_free(self.idmesh)



# ------------------------------------------------------------------------
def ElementParentFind(int                                           dnelt,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] elt_def_idx,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] elt_def,
                      int                                           dnelt_to_find,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] elt_to_find_def_idx,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] elt_to_find_def,
                      MPI.Comm    comm,      
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] parent):
    """
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    
    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    PDM_elt_parent_find(dnelt,
                        <int *>         elt_def_idx.data,
                        <PDM_g_num_t *> elt_def.data,
                        dnelt_to_find,
                        <int *>         elt_to_find_def_idx.data,
                        <PDM_g_num_t *> elt_to_find_def.data,
                        PDMC,     
                        <PDM_g_num_t *> parent.data)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::


# ------------------------------------------------------------------------
def ComputeDistributionFromDelmt(int         dnelt,
                                 MPI.Comm    comm, 
                                 int         offset=0): 
    """
    """
    # ************************************************************************
    # > Declaration
    cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='c'] elt_distrib = NPY.empty( comm.Get_size() + 1, dtype=npy_pdm_gnum_dtype, order='C')
    # ************************************************************************

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    
    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    PDM_distrib_compute(dnelt,
                        <PDM_g_num_t *> elt_distrib.data,
                        offset,
                        PDMC)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    
    return elt_distrib

# ------------------------------------------------------------------------

    
