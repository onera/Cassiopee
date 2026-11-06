# - getNormalMap (array) -
import Geom as D
import Generator as G
import Converter as C

# 2D structured
a = D.sphere((0,0,0), 1, 50)
n = G.getNormalMap(a)
n = C.center2Node(n);n = C.addVars([a, n])
C.convertArrays2File([n], 'out1.plt')

# 2D unstructured
a = D.sphere((0,0,0), 1, 50)
a = C.convertArray2Tetra(a)
n = G.getNormalMap(a)
n = C.center2Node(n);n = C.addVars([a, n])
C.convertArrays2File([n], 'out2.plt')
