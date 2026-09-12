# - close (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Transform.PyTree as T
import KCore.test as test

# BE NODE
a = G.cartHexa((0,0,0), (1,1,1), (5, 5, 5))
G._close(a)
test.testT(a, 1)

# ME NODE, no overlap
a = G.cartHexa((0,0,0), (1,1,1), (5, 5, 5))
b = G.cartHexa((4.5,0,0), (1,1,1), (5, 5, 5))
C._initVars(a, 'Density', 1.)
C._initVars(b, 'Density', 2.)
a = C.convertArray2Node(a)
b = C.convertArray2Node(b)
a = T.join(a, b)
C.convertPyTree2File(a, "out.cgns")
G._close(a)
test.testT(a, 2)

# ME NODE, overlap
a = G.cartHexa((0,0,0), (1,1,1), (5, 5, 5))
b = G.cartHexa((3,0,0), (1,1,1), (5, 5, 5))
C._initVars(a, 'Density', 1.)
C._initVars(b, 'Density', 2.)
a = C.convertArray2Node(a)
b = C.convertArray2Node(b)
a = T.join(a, b)
C.convertPyTree2File(a, "out.cgns")
G._close(a)
test.testT(a, 2)
