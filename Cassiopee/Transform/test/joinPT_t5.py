# - join (pyTree) -
import Transform.PyTree as T
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test
import numpy

# -- Join a list of BE
# 2D identical BE (quad-quad-quad) without fields at centers
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.), (5,10,1))
b = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.), (5,10,1))
c = G.cartHexa((0.8,0.,0.), (0.1,0.1,0.), (5,10,1))
C._initVars(a, 'F', 1.); C._initVars(b, 'F', 2.); C._initVars(c, 'F', 3.)
a = T.join([a, b, c]); t = C.newPyTree(["Base", a])
test.testT(t, 1)

# 2D BE (tri-quad-tri) without fields at centers
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.), (5,10,1))
b = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.), (5,10,1))
c = G.cartTetra((0.8,0.,0.), (0.1,0.1,0.), (5,10,1))
C._initVars(a, 'F', 1.); C._initVars(b, 'F', 2.); C._initVars(c, 'F', 3.)
a = T.join([a, b, c]); t = C.newPyTree(["Base", a])
test.testT(t, 2)

# 3D BE (pyra-hexa-penta-tetra) without fields at centers
a = G.cartPyra((0.,0.,0.), (0.1,0.1,0.1), (5,10,4))
b = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.1), (5,10,4))
c = G.cartPenta((0.8,0.,0.), (0.1,0.1,0.1), (5,10,4))
d = G.cartTetra((1.3,0.,0.), (0.1,0.1,0.1), (5,10,4))
C._initVars(a, 'F', 1.); C._initVars(b, 'F', 2.); C._initVars(c, 'F', 3.); C._initVars(d, 'F', 4.)
a = T.join([a, b, c, d]); t = C.newPyTree(["Base", a])
test.testT(t, 3)

# 3D BE (pyra-hexa-penta-hexa) with fields at centers
a = G.cartPyra((0.,0.,0.), (0.1,0.1,0.1), (5,10,4))
b = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.1), (5,10,4))
c = G.cartPenta((0.8,0.,0.), (0.1,0.1,0.1), (5,10,4))
d = G.cartHexa((1.2,0.,0.), (0.1,0.1,0.1), (5,10,4))
C._initVars(a, 'F', 1.); C._initVars(b, 'F', 2.); C._initVars(c, 'F', 3.); C._initVars(d, 'F', 4.)
C._initVars(a, 'centers:G', 4.); C._initVars(b, 'centers:G', 3.); C._initVars(c, 'centers:G', 2.); C._initVars(d, 'centers:G', 1.)
a = T.join([a, b, c, d]); t = C.newPyTree(["Base", a])
test.testT(t, 4)

# -- Join a list of NGON, api 3
# without fields at centers
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (5,10,4), api=3)
b = G.cartNGon((0.4,0.,0.), (0.1,0.1,0.1), (5,10,4), api=3)
c = G.cartNGon((0.8,0.,0.), (0.1,0.1,0.1), (5,10,4), api=3)
C._initVars(a, 'F', 1.); C._initVars(b, 'F', 2.); C._initVars(c, 'F', 3.)
a = T.join([a, b, c]); t = C.newPyTree(["Base", a])
test.testT(t, 5)

# with fields at centers
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (5,10,4), api=3)
b = G.cartNGon((0.4,0.,0.), (0.1,0.1,0.1), (5,10,4), api=3)
c = G.cartNGon((0.8,0.,0.), (0.1,0.1,0.1), (5,10,4), api=3)
C._initVars(a, 'F', 1.); C._initVars(b, 'F', 2.); C._initVars(c, 'F', 3.)
C._initVars(a, 'centers:G', 1.); C._initVars(b, 'centers:G', 2.); C._initVars(c, 'centers:G', 3.)
a = T.join([a, b, c]); t = C.newPyTree(["Base", a])
test.testT(t, 6)

# -- Join a list containing at least an empty connectivity
def createEmptyZone(etype):
    z = Internal.newZone(name="empty", zsize=[[0,0]], ztype="Unstructured")
    gc = Internal.newGridCoordinates(parent=z)
    Internal.newDataArray('CoordinateX', value=numpy.empty(0), parent=gc)
    Internal.newDataArray('CoordinateY', value=numpy.empty(0), parent=gc)
    Internal.newDataArray('CoordinateZ', value=numpy.empty(0), parent=gc)
    Internal.newElements(name="GridElements", etype=etype,
                         econnectivity=numpy.empty(0), erange=[1, 0],
                         eboundary=0, parent=z)
    return z

# Empty QUAD in second position ignored
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,10,1))
b = createEmptyZone(etype=7)
a = T.join([a, b]); t = C.newPyTree(["Base", a])
test.testT(t, 7)

# Empty HEXA in first position ignored
a = createEmptyZone(etype=17)
b = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,10,4))
a = T.join([a, b]); t = C.newPyTree(["Base", a])
test.testT(t, 8)

# Empty NODE with a valid HEXA
a = createEmptyZone(etype=2)
b = createEmptyZone(etype=2)
c = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,10,1))
d = createEmptyZone(etype=2)
a = T.join([a, b, c, d]); t = C.newPyTree(["Base", a])
test.testT(t, 9)

# Empty NODE with a valid NGON
a = createEmptyZone(etype=2)
b = createEmptyZone(etype=2)
c = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (5,10,1), api=3)
d = createEmptyZone(etype=2)
a = T.join([a, b, c, d]); t = C.newPyTree(["Base", a])
test.testT(t, 10)

# Two empty QUAD connectivities -> empty QUAD
a = createEmptyZone(etype=7)
b = createEmptyZone(etype=7)
a = T.join([a, b]); t = C.newPyTree(["Base", a])
test.testT(t, 11)

# Two empty NODE connectivities -> empty NODE
a = createEmptyZone(etype=2)
b = createEmptyZone(etype=2)
a = T.join([a, b]); t = C.newPyTree(["Base", a])
test.testT(t, 12)

# Two empty connectivities -> empty NODE
a = createEmptyZone(etype=5)
b = createEmptyZone(etype=7)
a = T.join([a, b]); t = C.newPyTree(["Base", a])
test.testT(t, 13)

# Two empty NGON connectivities -> empty NODE
a = createEmptyZone(etype=22)
b = createEmptyZone(etype=22)
a = T.join([a, b]); t = C.newPyTree(["Base", a])
test.testT(t, 14)
