# - zipper -
import Post as P
import Converter as C
import Transform as T
#
# extraction seulement des blocs aube + moyeu + conge
#
a = C.convertFile2Arrays("conge_complet.plt", "bin_tp")
arrays = [a[2], a[12],a[14],a[15]]
C.convertArrays2File(arrays, "out.plt", "bin_tp")

# Met cellN = 1 sur la fontiere (points normaux car on a supprime
# le domaine recouvert)
b = arrays[3]
s = T.subzone(b, (1,1,1), (b[2],1,b[4]))
s = C.initVars(s, 'cellN', 1)
b = T.patch(s, b, (1,1,1))
arrays[3] = b

# zipping
array = P.zipper(arrays, [])
C.convertArrays2File([array], "out2.plt","bin_tp")
