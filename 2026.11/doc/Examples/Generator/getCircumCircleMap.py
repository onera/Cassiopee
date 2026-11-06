# - getCircumCircleMap (array) -
import Geom as D
import Generator as G
import Converter as C

a = D.sphere((0,0,0), 1, 50)
a = C.convertArray2Tetra(a)
n = G.getCircumCircleMap(a)
n = C.center2Node(n); n = C.addVars([a, n])
C.convertArrays2File([n], "out.plt")
