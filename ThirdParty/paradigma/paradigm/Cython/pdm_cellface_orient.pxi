
cdef extern from "pdm_cellface_orient.h":

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    void PDM_cellface_orient(int           nCell,
                             int           nFace,
                             int           nVtx,
                             double*       coords,
                             int*          cellFaceIdx,
                             int*          cellFace,
                             int*          faceCell,
                             int*          faceVtxIdx,
                             int*          faceVtx)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


def PyPDM_cellface_orient(nCell, nFace, nVtx,
                          NPY.ndarray[NPY.double_t   , mode='c', ndim=1] VtxCoord not None,
                          NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] cellFaceIdx not None,
                          NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] cellFace not None,
                          NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] faceCell not None,
                          NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] faceVtxIdx not None,
                          NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] faceVtx not None
                          ):
  """
  Wrapping of PDM_cellface_orient
  """
  # print '*'*100
  # print 'nCell : ', nCell
  # print 'nFace : ', nFace
  # print 'nVtx  : ', nVtx
  # print 'VtxCoord  : ', VtxCoord
  # print 'cellFaceIdx  : ', cellFaceIdx
  # print 'cellFace  : ', cellFace
  # print 'faceVtxIdx  : ', faceVtxIdx
  # print 'faceVtx  : ', faceVtx
  PDM_cellface_orient(nCell,
                      nFace,
                      nVtx,
                      <double *> &VtxCoord[0],
                      <int *>    &cellFaceIdx[0],
                      <int *>    &cellFace[0],
                      <int *>    &faceCell[0],
                      <int *>    &faceVtxIdx[0],
                      <int *>    &faceVtx[0])


