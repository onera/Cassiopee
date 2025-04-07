# - zipper -
import Post as P
import Converter as C

# Maillage dauphin
arrays = C.convertFile2Arrays("dauphin_skin.plt","bin_tp")
array = P.zipper(arrays,['overlapTol',1.e-3])
C.convertArrays2File([array], "new.plt","bin_tp")
