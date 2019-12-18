
cdef extern from "pdm_part_to_block.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_part_to_block_t:
        pass

    ctypedef enum PDM_part_to_block_distrib_t:
        pass

    ctypedef enum PDM_part_to_block_post_t:
        pass

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_part_to_block_t *PDM_part_to_block_create(PDM_part_to_block_distrib_t   t_distrib,
                                                  PDM_part_to_block_post_t      t_post,
                                                  double                         partActiveNode,
                                                  PDM_g_num_t                 **gnum_elt,
                                                  double                       **weight,
                                                  int                          *n_elt,
                                                  int                           n_part,
                                                  PDM_MPI_Comm                  comm)

    int PDM_part_to_block_n_active_ranks_get(PDM_part_to_block_t *ptb)

    int *PDM_part_to_block_active_ranks_get(PDM_part_to_block_t *ptb)

    int PDM_part_to_block_is_active_rank(PDM_part_to_block_t *ptb)

    int PDM_part_to_block_n_elt_block_get(PDM_part_to_block_t *ptb)

    PDM_g_num_t *PDM_part_to_block_block_gnum_get(PDM_part_to_block_t *ptb)

    int PDM_part_to_block_exch(PDM_part_to_block_t       *ptb,
                               size_t                     s_data,
                               PDM_stride_t               t_stride,
                               int                        cst_stride,
                               int                      **part_stride,
                               void                     **part_data,
                               int                      **block_stride,
                               void                     **block_data)

    PDM_part_to_block_t *PDM_part_to_block_free(PDM_part_to_block_t *ptb)

    PDM_g_num_t *PDM_part_to_block_distrib_index_get(PDM_part_to_block_t *ptb)

    PDM_l_num_t *PDM_part_to_block_destination_get(PDM_part_to_block_t *ptb)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# ------------------------------------------------------------------
cdef class PartToBlock:
    """
       PartToBlock: Interface for block_to_part.c
    """
    # ************************************************************************
    # > Class attributes
    cdef PDM_part_to_block_t *PTB
    cdef int                  partN
    cdef int                  Size
    cdef int                  Rank

    cdef int                 *NbElmts
    cdef PDM_g_num_t        **LNToGN
    cdef double              **weight

    cdef PDM_part_to_block_distrib_t t_distrib
    cdef PDM_part_to_block_post_t    t_post
    cdef PDM_stride_t                t_stride

# ************************************************************************
    # ------------------------------------------------------------------------
    def __cinit__(self, MPI.Comm comm, list pLNToGN, list pWeight, int partN,
                        PDM_part_to_block_distrib_t t_distrib = <PDM_part_to_block_distrib_t> (0),
                        PDM_part_to_block_post_t    t_post    = <PDM_part_to_block_post_t   > (0),
                        PDM_stride_t                t_stride  = <PDM_stride_t   > (0),
                        double partActiveNode = 1.):
        """
        TODOUX
        """
        # ************************************************************************
        # > Declaration
        cdef int      nElts
        cdef int      idx
        # > Numpy array
        cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='fortran'] partLNToGN
        cdef NPY.ndarray[double, ndim=1, mode='fortran'] partWeight
        # ************************************************************************

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Some verification
        assert(len(pLNToGN) == partN)
        if (pWeight is not None):
          assert(len(pWeight) == partN)

        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Store partN and parameter
        self.partN     = partN
        self.t_distrib = t_distrib
        self.t_post    = t_post
        self.t_stride  = t_stride

        self.Rank = comm.Get_rank()
        self.Size = comm.Get_size()
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Convert mpi4py -> PDM_MPI
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi
        cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Allocate
        self.LNToGN   = <PDM_g_num_t **> malloc(sizeof(PDM_g_num_t *) * self.partN )
        self.weight   = NULL

        if (pWeight is not None):
          self.weight   = <double **> malloc(sizeof(double *) * self.partN )


        self.NbElmts  = <int * > malloc(sizeof(int  ) * self.partN )
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Prepare
        for idx, partLNToGN in enumerate(pLNToGN):

          # ------------------------------------------------
          # > Get shape of array
          nElts = partLNToGN.shape[0]
          self.NbElmts[idx] = <int> nElts
          # ------------------------------------------------

          # ------------------------------------------------
          # > Assign array
          self.LNToGN[idx] = <PDM_g_num_t *> partLNToGN.data
          # ------------------------------------------------
          # print "nElts, partLNToGN",nElts, partLNToGN

        idx = 0
        if (pWeight is not None):
          for idx, partWeight in enumerate(pWeight):
            self.weight[idx] = <double *> partWeight.data


        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Create
        self.PTB = PDM_part_to_block_create(t_distrib,
                                            t_post,
                                            partActiveNode,
                                            self.LNToGN,
                                            self.weight,
                                            self.NbElmts,
                                            self.partN,
                                            PDMC)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def PartToBlock_Exchange(self, dict dField, dict pField, list pStrid = None):
        """
           TODOUX : 1) Exchange of variables types array
                    2) Assertion of type and accross MPI of the same field
        """
        # ************************************************************************
        # > Declaration
        cdef NPY.ndarray   dArray
        cdef NPY.ndarray   pArray
        cdef NPY.ndarray   partStrid
        cdef int           idx
        cdef int           blk_size
        cdef NPY.npy_intp *ArrayDim

        # > For PDM
        cdef size_t   s_data
        cdef int      strideOne
        cdef int    **part_stride
        cdef int     *block_stride
        cdef void    *block_data
        cdef void   **part_data
        cdef int      ndim
        cdef int      npyflags=-1;
        cdef NPY.ndarray tmpData
        # ************************************************************************

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Allocate
        part_data = <void **> malloc(self.partN * sizeof(void **))
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Prepare stride
        if(self.t_stride == 0): # Cst Stride
           strideOne     = 1
           part_stride   = NULL
        else:
           strideOne     = 0
           part_stride   = <int **> malloc(self.partN * sizeof(int *))
           assert(pStrid is not None)
           assert(len(pStrid) == self.partN)
           # for idx in xrange(self.partN):
           for idx, partStrid in enumerate(pStrid):
              part_stride[idx] = <int *> partStrid.data
        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        ndim = -1

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Loop on all field to build
        for field, partList in pField.iteritems():

          # print field, partList

          # ::::::::::::::::::::::::::::::::::::::::::::::::::
          # > Prepare part_data
          for idx, pArray in enumerate(partList):
            # ------------------------------------------------
            # > Get flow solution
            if(pArray.ndim == 2):
              assert(pArray.shape[1] == self.NbElmts[idx])
              ndim = 2
            else:
              if(self.t_stride == <PDM_stride_t   >(0)):
                assert(pArray.shape[0] == self.NbElmts[idx])
              ndim = 1
            part_data[idx] = <void *> pArray.data
            # ------------------------------------------------

            # ------------------------------------------------
            # > Fill s_data - How to check if array is different ?
            s_data     = pArray.dtype.itemsize
            dtype_data = pArray.dtype.num
            dtypep     = pArray.dtype
            # ------------------------------------------------

          # ::::::::::::::::::::::::::::::::::::::::::::::::::
          # > Prepare block_data
          block_stride  = NULL
          block_data    = NULL
          # ::::::::::::::::::::::::::::::::::::::::::::::::::

          # ::::::::::::::::::::::::::::::::::::::::::::::::::
          # > Exchange
          # print "PDM_part_to_block_exch "
          c_size = PDM_part_to_block_exch(self.PTB,
                                          s_data,
                                          self.t_stride,
                                          strideOne,
                                          part_stride,
                                          part_data,
                                          &block_stride,
                                          <void **> &block_data)
          # print "PDM_part_to_block_exch end "
          # ::::::::::::::::::::::::::::::::::::::::::::::::::

          # ::::::::::::::::::::::::::::::::::::::::::::::::::
          blk_size = PDM_part_to_block_n_elt_block_get(self.PTB);
          # ::::::::::::::::::::::::::::::::::::::::::::::::::

          # ::::::::::::::::::::::::::::::::::::::::::::::::::
          # > Put in dict
          dim           = <NPY.npy_intp> c_size
          if(c_size == 0):
            # dField[field] = None # Attention faire une caspule vide serait mieux non ?
            dField[field] = NPY.empty((0), dtype=dtypep)
            # print 'Attention in PDM_part_to_block'
          else:
            if(ndim == 2):
              ArrayDim    = <NPY.npy_intp *> malloc(2 * sizeof(NPY.npy_intp *))
              ArrayDim[0] = <NPY.npy_intp> 1
              ArrayDim[1] = <NPY.npy_intp> c_size

              # > Put in dField
              tmpData = NPY.PyArray_SimpleNewFromData(ndim, ArrayDim, dtype_data, <void *> block_data)
              PyArray_ENABLEFLAGS(tmpData, NPY.NPY_OWNDATA);
              dField[field] = tmpData
              # > Free
              free(ArrayDim)
            else:
              # dField[field] = NPY.PyArray_SimpleNewFromData(1, &dim, dtype_data,
              #                                               <void *> block_data)
              tmpData = NPY.PyArray_SimpleNewFromData(1, &dim, dtype_data,
                                                         <void *> block_data)
              PyArray_ENABLEFLAGS(tmpData, NPY.NPY_OWNDATA);
              dField[field] = tmpData
              # print(dField[field].flags)
          # ::::::::::::::::::::::::::::::::::::::::::::::::::

          # ::::::::::::::::::::::::::::::::::::::::::::::::::
          # > Stride management
          if(self.t_stride == 1 ):
            dimStri = <NPY.npy_intp> blk_size
            # dField[field+'#Stride'] = NPY.PyArray_SimpleNewFromData(1, &dimStri, NPY.NPY_INT32, <void *> block_stride)
            tmpData = NPY.PyArray_SimpleNewFromData(1, &dimStri, NPY.NPY_INT32, <void *> block_stride)
            PyArray_ENABLEFLAGS(tmpData, NPY.NPY_OWNDATA);
            dField[field+'#Stride'] = tmpData
          # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Deallocate
        free(part_data)
        if(self.t_stride != 0): # Var Stride
          free(part_stride)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def getDistributionCopy(self):
      """
         Return a copy of the distrisbution array compute in library
         Copy because remove of PTB object can made a core ...
      """
      # ************************************************************************
      # > Declaration
      cdef PDM_g_num_t* Distrib
      # ************************************************************************

      # ::::::::::::::::::::::::::::::::::::::::::::::::::
      # > Get
      Distrib = PDM_part_to_block_distrib_index_get(self.PTB)
      # ::::::::::::::::::::::::::::::::::::::::::::::::::

      # ::::::::::::::::::::::::::::::::::::::::::::::::::
      dim        = <NPY.npy_intp> (self.Size+1)
      DistribNPY = NPY.PyArray_SimpleNewFromData(1, &dim, PDM_G_NUM_NPY_INT, <void *> Distrib)
      # ::::::::::::::::::::::::::::::::::::::::::::::::::

      return NPY.copy(DistribNPY)

    # ------------------------------------------------------------------------
    def getBeginNbEntryAndGlob(self):
      """
         Return a copy of the distrisbution array compute in library
         Copy because remove of PTB object can made a core ...
      """
      # ************************************************************************
      # > Declaration
      cdef PDM_g_num_t* Distrib
      # ************************************************************************

      # ::::::::::::::::::::::::::::::::::::::::::::::::::
      # > Get
      Distrib = PDM_part_to_block_distrib_index_get(self.PTB)
      # ::::::::::::::::::::::::::::::::::::::::::::::::::

      # ::::::::::::::::::::::::::::::::::::::::::::::::::
      dim        = <NPY.npy_intp> (self.Size+1)
      DistribNPY = NPY.PyArray_SimpleNewFromData(1, &dim, PDM_G_NUM_NPY_INT, <void *> Distrib)
      # ::::::::::::::::::::::::::::::::::::::::::::::::::

      Beg = DistribNPY[self.Rank]
      NbE = DistribNPY[self.Rank+1]-DistribNPY[self.Rank]
      GlB = DistribNPY[self.Size]

      return (Beg, NbE, GlB)

    # ------------------------------------------------------------------------
    def __dealloc__(self):
      """
         Use the free method of PDM Lib
      """
      # ************************************************************************
      # > Declaration
      cdef PDM_part_to_block_t *a
      # ************************************************************************

      # > Free Ppart Structure
      a = PDM_part_to_block_free(self.PTB)

      # > Free allocated array
      free(self.LNToGN)
      if (self.weight != NULL):
        free(self.weight)
      free(self.NbElmts)

