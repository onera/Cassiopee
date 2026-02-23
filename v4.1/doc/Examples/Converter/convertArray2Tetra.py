# - convertArray2Tetra (array) -
import Converter as C
import Generator as G

# 2D: triangles
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
b = C.convertArray2Tetra(a, split="withBarycenters")
C.convertArrays2File([b], "out.plt")
