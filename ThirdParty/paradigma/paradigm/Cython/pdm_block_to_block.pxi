
cdef extern from "pdm_block_to_block.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure 
    ctypedef struct PDM_block_to_block_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function 
    PDM_block_to_block_t *PDM_block_to_block_create(PDM_g_num_t   *blockDistribIniIdx,
                                                    PDM_g_num_t   *blockDistribEndIdx,
                                                    PDM_MPI_Comm   comm)

    int PDM_block_to_block_exch(PDM_block_to_block_t  *btb,
                                size_t                s_data,
                                PDM_stride_t          t_stride,
                                int                  *block_stride_ini,
                                void                 *block_data_ini,
                                int                  *block_stride_end,
                                void                **block_data_end)
    
    int PDM_block_to_block_exch_int(PDM_block_to_block_t  *btb,
                                    size_t                s_data,
                                    PDM_stride_t          t_stride,
                                    int                  *block_stride_ini,
                                    int                  *block_data_ini,
                                    int                  *block_stride_end,
                                    int                 **block_data_end)

    PDM_block_to_block_t *PDM_block_to_block_free(PDM_block_to_block_t *btb)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# ------------------------------------------------------------------
cdef class BlockToBlock:
    """
       BlockToPart: Interface for block_to_block.c
    """
    # ************************************************************************
    # > Class attributes
    cdef PDM_block_to_block_t *BTB
    # cdef int  *BlkDistribIdx
    # ************************************************************************
    # ------------------------------------------------------------------------
    def __cinit__(self, NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='c'] DistribIni,
                        NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='c'] DistribEnd,
                        MPI.Comm comm):
        """
        Constructor of BlockToBlock object : Python wrapping of PDM library (E. Quémerais)

            :param comm:        MPI Communicator (Caution MPI Comm is a mpi4py object )
            :param DistribIni:  Distribution of distribute array in initiale frame (Size = nRank+1)
            :param DistribEnd:  Distribution of distribute array in finale frame (Size = nRank+1)

        """
        # ************************************************************************
        # > Declaration
        # ************************************************************************
        
        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Some verification
        assert(DistribIni.shape[0] == DistribEnd.shape[0])
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Mpi size
        Size = comm.Get_size()
        assert(DistribIni.shape[0] == Size+1)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Convert mpi4py -> PDM_MPI
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi
        cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Create
        # self.btb = PDM_block_to_block_create(self.BlkDistribIdx, self.LNToGN,
        self.BTB = PDM_block_to_block_create(<PDM_g_num_t *> DistribIni.data, 
                                             <PDM_g_num_t *> DistribEnd.data, 
                                             PDMC)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------
    def BlockToBlock_Exchange(self, dict         dFieldIni, 
                                    dict         dFieldEnd,
                                    PDM_stride_t t_stride = <PDM_stride_t> 0):
        """
           TODOUX : 1) Exchange of variables types array 
                    2) Assertion of type and accross MPI of the same field
                    3) Stride variable !!
        """
        # ************************************************************************
        # > Declaration
        cdef NPY.ndarray dArrayIni
        cdef NPY.ndarray dArrayEnd

        # > Specific to PDM
        cdef size_t      s_data
        cdef int        *block_stride_ini
        cdef int        *block_stride_end
        cdef void       *block_data_ini
        cdef void       *block_data_end
        # ************************************************************************

        # > Understand how to interface
        block_stride_ini = <int *> malloc( 1 * sizeof(int *))
        block_stride_end = <int *> malloc( 1 * sizeof(int *))

        # > Init
        block_stride_ini[0] = 1 # No entrelacing data
        block_stride_end[0] = 1 # No entrelacing data

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Loop over all field and fill buffer
        #   En MPI il y a un risque que les dictionnaires ne soit pas les meme ?
        #   Ordered Dict ??
        # mpi.gather puis si vrai 0 et assert sur la somme ?
        for field in dFieldIni.keys():
          # ::::::::::::::::::::::::::::::::::::::::::::::::::
          # > Get field and assign ptr
          dArrayIni      = dFieldIni[field]
          block_data_ini = <void *> dArrayIni.data
          block_data_end = NULL
          s_data         = dArrayIni.dtype.itemsize
          dtype_data     = dArrayIni.dtype.num
          # ::::::::::::::::::::::::::::::::::::::::::::::::::

          # ::::::::::::::::::::::::::::::::::::::::::::::::::
          # > Compute
          # PDM_block_to_block_exch(self.BTB, 
          #                         s_data, 
          #                         t_stride,
          #                         block_stride_ini,
          #                         block_data_ini,
          #                         block_stride_end,
          #                         block_data_end)
          c_size = PDM_block_to_block_exch_int(self.BTB, 
                                               s_data, 
                                               t_stride,
                                               block_stride_ini,
                                               <int *> block_data_ini,
                                               block_stride_end,
                                               <int **> &block_data_end)
          
          dim = <NPY.npy_intp> c_size
          # print "c_size : ", c_size
          if(c_size == 0):
            dFieldEnd[field] = None
          else:
            dFieldEnd[field] = NPY.PyArray_SimpleNewFromData(1, &dim, dtype_data, <void *> block_data_end)
            # print dFieldEnd[field]
          # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------
    def __dealloc__(self):
      """ 
         Deallocate all the array 
      """
      # ************************************************************************
      # > Declaration
      cdef PDM_block_to_block_t *a
      # ************************************************************************
      # > Free PDM Structure
      a = PDM_block_to_block_free(self.BTB)

