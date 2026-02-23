# - splitSharpEdges (array) -
import Converter as C
import Transform as T
import Geom as D
import Generator as G

a = D.text3D("A"); a = G.close(a, 1.e-4)
B = T.splitSharpEdges(a, 89.)
C.convertArrays2File(B, 'out.plt')
