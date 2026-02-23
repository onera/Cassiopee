# - splitSharpEdges (pyTree) -
import Converter.PyTree as C
import Transform.PyTree as T
import Geom.PyTree as D
import Generator.PyTree as G
import KCore.test as test

a = D.text3D("A", font='text1'); a = G.close(a, 1.e-3)
B = T.splitSharpEdges(a, 89.)
t = C.newPyTree(['Base',2]); t[2][1][2] += B
test.testT(t, 1)
