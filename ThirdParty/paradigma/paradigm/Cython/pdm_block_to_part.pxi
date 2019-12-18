
cdef extern from "pdm_block_to_part.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_block_to_part_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_block_to_part_t *PDM_block_to_part_create(PDM_g_num_t   *blockDistribIdx,
                                                  PDM_g_num_t  **gnum_elt,
                                                  int           *n_elt,
                                                  int            n_part,
                                                  PDM_MPI_Comm   comm)

    void PDM_block_to_part_exch(PDM_block_to_part_t  *btp,
                                size_t                s_data,
                                PDM_stride_t          t_stride,
                                int                  *block_stride,
                                void                 *block_data,
                                int                 **part_stride,
                                void                **part_data)

    void PDM_block_to_part_exch2(PDM_block_to_part_t  *btp,
                                size_t                s_data,
                                PDM_stride_t          t_stride,
                                int                  *block_stride,
                                void                 *block_data,
                                int                 ***part_stride,
                                void                ***part_data)

    PDM_block_to_part_t *PDM_block_to_part_free(PDM_block_to_part_t *btp)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# ------------------------------------------------------------------
cdef class BlockToPart:
    """
       BlockToPart: Interface for block_to_part.c
    """
    # ************************************************************************
    # > Class attributes
    cdef PDM_block_to_part_t *BTP
    cdef int                  partN
    # > Assertion and debug (To be removed )
    cdef int          *NbElmts
    cdef PDM_g_num_t **LNToGN
    # cdef int  *BlkDistribIdx
    # ************************************************************************
    # ------------------------------------------------------------------------
    def __cinit__(self, NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='c'] Distrib,
                        MPI.Comm comm,
                        list     pLNToGN,
                        int      partN):
        """
        Constructor of BlockToPart object : Python wrapping of PDM library (E. Quémerais)

            :param comm:     MPI Communicator (Caution MPI Comm is a mpi4py object )
            :param Distrib:  Distribution of distribute array (Size = nRank+1)
            :param pLNToGN:  Part list containaing numpy on LNToGN for each partition (len = partN)
            :param partN:    Number of partition

        """
        # ************************************************************************
        # > Declaration
        cdef int      nElts
        cdef int      idx
        # > Numpy array
        cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='fortran'] partLNToGN
        # ************************************************************************
        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Some verification
        assert(len(pLNToGN) == partN)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Store partN
        self.partN = partN
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Mpi size
        Size = comm.Get_size()
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Convert mpi4py -> PDM_MPI
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi
        cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Allocate
        self.LNToGN        = <PDM_g_num_t **> malloc(sizeof(PDM_g_num_t **) * partN )
        self.NbElmts       = <int *         > malloc(sizeof(int *         ) * partN )
        # self.BlkDistribIdx = <int * > malloc(sizeof(int * ) * Size+1)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Assign
        # self.BlkDistribIdx = <int *>  Distrib.data
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
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Create
        # self.BTP = PDM_block_to_part_create(self.BlkDistribIdx, self.LNToGN,
        self.BTP = PDM_block_to_part_create(<PDM_g_num_t *> Distrib.data,
                                            self.LNToGN,
                                            self.NbElmts,
                                            self.partN,
                                            PDMC)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------
    def BlockToPart_Exchange(self, dict         dField,
                                   dict         pField,
                                   PDM_stride_t t_stride = <PDM_stride_t> 0,
                                   NPY.ndarray[NPY.int32_t, ndim=1, mode='c'] BlkStride = None):
        """
           TODOUX : 1) Exchange of variables types array
                    2) Assertion of type and accross MPI of the same field
        """
        # ************************************************************************
        # > Declaration
        cdef NPY.ndarray dArray
        cdef NPY.ndarray pArray

        # > Specific to PDM
        cdef size_t      s_data
        cdef int        *block_stride
        cdef int       **part_stride
        cdef void      **part_data
        cdef void       *block_data
        # ************************************************************************

        # > Understand how to interface
        if(BlkStride is None):
          block_stride = <int *> malloc( 1 * sizeof(int *))
          part_stride  = NULL
          part_data    = <void **> malloc(self.partN * sizeof(void **))

          # > Init
          block_stride[0] = 1 # No entrelacing data
        else:
          block_stride = <int *> BlkStride.data
          part_stride  = <int **>  malloc(self.partN * sizeof(void **))
          part_data    = <void **> malloc(self.partN * sizeof(void **))
          # part_stride  = NULL
          # part_data    = NULL

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Loop over all field and fill buffer
        #   En MPI il y a un risque que les dictionnaires ne soit pas les meme ?
        #   Ordered Dict ??
        # mpi.gather puis si vrai 0 et assert sur la somme ?
        for field in dField.keys():
          # ::::::::::::::::::::::::::::::::::::::::::::::::::
          # > Get field and assign ptr
          dArray     = dField[field]
          block_data = <void *> dArray.data
          s_data     = dArray.dtype.itemsize

          # ::::::::::::::::::::::::::::::::::::::::::::::::::
          # > Get part list and Loop over all part
          partList = pField[field]
          for idx, pArray in enumerate(partList):
            # if(pArray.ndim == 2):
            #   assert(pArray.shape[1] == self.NbElmts[idx])
            # else:
            #   assert(pArray.shape[0] == self.NbElmts[idx])
            part_data[idx] = <void *> pArray.data

          if(BlkStride is not None):
            partStrideList = pField[field+'#PDM_Stride']
            for idx, pArray in enumerate(partStrideList):
              part_stride[idx] = <int *> pArray.data

          # ::::::::::::::::::::::::::::::::::::::::::::::::::
          # > Compute
          PDM_block_to_part_exch(self.BTP,
                                 s_data,
                                 t_stride,
                                 block_stride,
                                 block_data,
                                 part_stride,
                                 part_data)
          # > Verbose
          # print field
          # for idx, pArray in enumerate(partList):
          #   print pArray

        # > Deallocate - On dessaloue uniquement les pointeurs de partition
        #   Le coté pointeur de pointeur est assuré par la list python
        if(BlkStride is None):
          free(block_stride)
          free(part_data)
        else:
          free(block_stride)
          free(part_data)
          free(part_stride)


    # ------------------------------------------------------------------
    def BlockToPart_Exchange2(self, dict         dField,
                                    dict         pField,
                                    PDM_stride_t t_stride = <PDM_stride_t> 0,
                                    NPY.ndarray[NPY.int32_t, ndim=1, mode='c'] BlkStride = None):
        """
           TODOUX : 1) Exchange of variables types array
                    2) Assertion of type and accross MPI of the same field
        """
        # ************************************************************************
        # > Declaration
        cdef NPY.ndarray dArray
        cdef NPY.ndarray pArray

        # > Specific to PDM
        cdef int         ii, sizeData
        cdef size_t      s_data
        cdef int        *block_stride
        cdef int       **part_stride
        cdef void      **part_data
        cdef void       *block_data
        cdef NPY.ndarray tmpData
        cdef NPY.ndarray tmpStride
        # ************************************************************************

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        for field in dField.keys():
          # ::::::::::::::::::::::::::::::::::::::::::::::::::
          # > Get field and assign ptr
          dArray     = dField[field]
          block_data = <void *> dArray.data
          s_data     = dArray.dtype.itemsize

          # > Understand how to interface
          if(BlkStride is None):
            block_stride = <int *> malloc( 1 * sizeof(int *))
            part_stride  = NULL
            part_data    = NULL
            # > Init
            block_stride[0] = 1 # No entrelacing data
          else:
            block_stride = <int *> BlkStride.data
            part_stride  = NULL
            part_data    = NULL
          # ::::::::::::::::::::::::::::::::::::::::::::::::::

          # ::::::::::::::::::::::::::::::::::::::::::::::::::
          # > Compute
          PDM_block_to_part_exch2(self.BTP,
                                  s_data,
                                  t_stride,
                                  block_stride,
                                  block_data,
                                  &part_stride,
                                  &part_data)
          # ::::::::::::::::::::::::::::::::::::::::::::::::::

          # ::::::::::::::::::::::::::::::::::::::::::::::::::
          if(BlkStride is not None):
            pField[field] = list()
            pField[field+"#PDM_Stride"] = list()
            for idx in xrange(self.partN):
              sizeData = 0
              for ii in xrange(self.NbElmts[idx]):
                sizeData += part_stride[idx][ii]

              # print("sizeData[{0}] = {1}".format(idx, sizeData))
              dimData = <NPY.npy_intp> sizeData
              # Memory leak block_data will be never desalocated ...
              tmpData = NPY.PyArray_SimpleNewFromData(1, &dimData, dArray.dtype.num, <void *> part_data[idx])
              PyArray_ENABLEFLAGS(tmpData, NPY.NPY_OWNDATA);
              pField[field].append(tmpData)

              dimStri = <NPY.npy_intp> self.NbElmts[idx]
              # print("dimStri : ", dimStri)
              tmpStride = NPY.PyArray_SimpleNewFromData(1, &dimStri, NPY.NPY_INT32, <void *> part_stride[idx])
              PyArray_ENABLEFLAGS(tmpStride, NPY.NPY_OWNDATA);
              pField[field+"#PDM_Stride"].append(tmpStride)
          # ::::::::::::::::::::::::::::::::::::::::::::::::::

          # ::::::::::::::::::::::::::::::::::::::::::::::::::
          free(part_data)
          free(part_stride)
          # ::::::::::::::::::::::::::::::::::::::::::::::::::
          # free(part_stride) # Les tableaux en eux meme sont detenu par le python

    # ------------------------------------------------------------------
    def __dealloc__(self):
      """
         Deallocate all the array
      """
      # ************************************************************************
      # > Declaration
      cdef PDM_block_to_part_t *a
      # ************************************************************************
      # > Free PDM Structure
      a = PDM_block_to_part_free(self.BTP)

      # > Free allocated array
      free(self.NbElmts)
      free(self.LNToGN)
      # free(self.BlkDistribIdx)s

