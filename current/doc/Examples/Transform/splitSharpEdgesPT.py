# - splitSharpEdges (pyTree) -
import Converter.PyTree as C
import Transform.PyTree as T
import Geom.PyTree as D
import Generator.PyTree as G

a = D.text3D("A"); a = G.close(a, 1.e-3)
B = T.splitSharpEdges(a, 89.)
C.convertPyTree2File(B, 'out.cgns')
