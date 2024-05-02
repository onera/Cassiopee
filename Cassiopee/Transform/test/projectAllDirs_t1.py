# - projectAllDirs (array) -
import Geom as D
import Converter as C
import Generator as G
import Transform as T
import KCore.test as test

# Structure
a = D.sphere((0,0,0), 1., 20)
a = C.initVars(a, 'F', 1)
b = G.cart((1.1,-0.1,-0.1),(0.1,0.1,0.1), (1,5,5))
n = G.getNormalMap(b); n = C.center2Node(n)
C._addVars([b,n])
b = C.initVars(b, 'F', 1)
c = T.projectAllDirs([b], [a], ['sx','sy','sz'])
test.testA(c, 1)

# Non structure
a = D.sphere((0,0,0), 1., 20)
a = C.initVars(a, 'F', 1)
b = G.cartTetra((1.1,-0.1,-0.1),(0.1,0.1,0.1), (1,5,5))
n = G.getNormalMap(b); n = C.center2Node(n)
C._addVars([b,n])
b = C.initVars(b, 'F', 1)
c = T.projectAllDirs([b], [a], ['sx','sy','sz'])
test.testA(c, 2)

a = C.convertArray2Tetra(a)
b = G.cartTetra((1.1,-0.1,-0.1),(0.1,0.1,0.1), (1,5,5))
n = G.getNormalMap(b); n = C.center2Node(n)
C._addVars([b,n])
b = C.initVars(b, 'F', 1)
c = T.projectAllDirs([b], [a], ['sx','sy','sz'])
test.testA(c, 3)
