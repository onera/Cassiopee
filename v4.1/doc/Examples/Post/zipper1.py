# - zipper -
import Post as P
import Converter as C

# Maillage spoiler
arrays = C.convertFile2Arrays("spoiler.plt", "bin_tp")
array = P.zipper(arrays, ['overlapTol',1.e-3])
C.convertArrays2File([array], "out.plt", "bin_tp")
