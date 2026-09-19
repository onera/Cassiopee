# - join (pyTree) -
import Transform.PyTree as T
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# --- Join two BE
# 2D identical BE (quad-quad) without fields at centers
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.), (5,10,1))
b = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.), (5,10,1))
C._initVars(a, 'F', 1.); C._initVars(b, 'F', 2.)
a = T.join(a, b); t = C.newPyTree(["Base", a])
test.testT(t, 1)

# 2D BE (tri-quad) without fields at centers
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.), (5,10,1))
b = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.), (5,10,1))
C._initVars(a, 'F', 1.); C._initVars(b, 'F', 2.)
a = T.join(a, b); t = C.newPyTree(["Base", a])
test.testT(t, 2)

# 3D BE (pyra-hexa) without fields at centers
a = G.cartPyra((0.,0.,0.), (0.1,0.1,0.1), (5,10,4))
b = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.1), (5,10,4))
C._initVars(a, 'F', 1.); C._initVars(b, 'F', 2.)
a = T.join(a, b); t = C.newPyTree(["Base", a])
test.testT(t, 3)

# 3D BE (pyra-hexa) with fields at centers
a = G.cartPyra((0.,0.,0.), (0.1,0.1,0.1), (5,10,4))
b = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.1), (5,10,4))
C._initVars(a, 'F', 1.); C._initVars(b, 'F', 2.)
C._initVars(a, 'centers:G', 3.); C._initVars(b, 'centers:G', 2.)
a = T.join(a, b); t = C.newPyTree(["Base", a])
test.testT(t, 4)

# --- Join two ME
# without fields at centers: tri-quad-tri
#                                 |
#                                tri
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (5,10,1))
b = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.2), (5,10,1))
a = C.mergeConnectivity([a, b])
c = G.cartTetra((0.8,0.,0.), (0.1,0.1,0.2), (5,10,1))
d = G.cartTetra((0.4,-0.9,0.), (0.1,0.1,0.2), (5,10,1))
c = C.mergeConnectivity([c, d])
C._initVars(a, 'F', 1.); C._initVars(c, 'F', 2.)
a = T.join(a, c); t = C.newPyTree(["Base", a])
test.testT(t, 5)

# with fields at centers: hexa - pyra
#                          |      |
#                        penta - hexa
a = G.cartHexa((0.,0.4,0.), (0.1,0.1,0.1), (5,5,5))
b = G.cartPyra((0.4,0.4,0.), (0.1,0.1,0.1), (5,5,5))
C._initVars(a, 'F', 1.); C._initVars(b, 'F', 2.)
C._initVars(a, 'centers:G', 1.); C._initVars(b, 'centers:G', 2.)
a = C.mergeConnectivity([a, b])
c = G.cartPenta((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
d = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.1), (5,5,5))
C._initVars(c, 'F', 3.); C._initVars(d, 'F', 4.)
C._initVars(c, 'centers:G', 3.); C._initVars(d, 'centers:G', 4.)
c = C.mergeConnectivity([c, d])
a = T.join(a, c); t = C.newPyTree(["Base", a])
#C.convertPyTree2File(t, "out.cgns"); exit()
#test.testT(t, 6)

# with fields at centers: hexa - pyra
#     (formula)            |      |
#                        penta - hexa
a = G.cartHexa((0.,0.4,0.), (0.1,0.1,0.1), (5,5,5))
b = G.cartPyra((0.4,0.4,0.), (0.1,0.1,0.1), (5,5,5))
C._initVars(a, 'F={CoordinateY}'); C._initVars(b, 'F={CoordinateY}')
C._initVars(a, 'centers:G={centers:CoordinateX}'); C._initVars(b, 'centers:G={centers:CoordinateX}')
a = C.mergeConnectivity([a, b])
c = G.cartPenta((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
d = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.1), (5,5,5))
C._initVars(c, 'F={CoordinateY}'); C._initVars(d, 'F={CoordinateY}')
C._initVars(c, 'centers:G={centers:CoordinateX}'); C._initVars(d, 'centers:G={centers:CoordinateX}')
c = C.mergeConnectivity([c, d])
a = T.join(a, c); t = C.newPyTree(["Base", a])
#C.convertPyTree2File(t, "out.cgns"); exit()
#test.testT(t, 7)

# --- Join two NGON, api 3
# without fields at centers
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (5,10,4), api=3)
b = G.cartNGon((0.4,0.,0.), (0.1,0.1,0.1), (5,10,4), api=3)
C._initVars(a, 'F', 1.); C._initVars(b, 'F', 2.)
a = T.join(a, b); t = C.newPyTree(["Base", a])
test.testT(t, 8)

# with fields at centers
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (5,10,4), api=3)
b = G.cartNGon((0.4,0.,0.), (0.1,0.1,0.1), (5,10,4), api=3)
C._initVars(a, 'F', 1.); C._initVars(b, 'F', 2.)
C._initVars(a, 'centers:G', 1.); C._initVars(b, 'centers:G', 2.)
a = T.join(a, b); t = C.newPyTree(["Base", a])
test.testT(t, 9)
