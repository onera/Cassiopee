# - join (pyTree) -
import Transform.PyTree as T
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# --- Join a list of BE
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
C._initVars(a, 'centers:G', 1.); C._initVars(b, 'centers:G', 2.); C._initVars(c, 'centers:G', 3.); C._initVars(d, 'centers:G', 4.)
a = T.join([a, b, c, d]); t = C.newPyTree(["Base", a])
#C.convertPyTree2File(t, "out.cgns"); exit()
#test.testT(t, 4)

# --- Join a list of NGON, api 3
# without fields at centers
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (5,10,4), api=3)
b = G.cartNGon((0.4,0.,0.), (0.1,0.1,0.1), (5,10,4), api=3)
c = G.cartNGon((0.8,0.,0.), (0.1,0.1,0.1), (5,10,4), api=3)
C._initVars(a, 'F', 1.); C._initVars(b, 'F', 2.); C._initVars(c, 'F', 3.)
a = T.join([a, b, c]); t = C.newPyTree(["Base", a])
#test.testT(t, 5)

# with fields at centers
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (5,10,4), api=3)
b = G.cartNGon((0.4,0.,0.), (0.1,0.1,0.1), (5,10,4), api=3)
c = G.cartNGon((0.8,0.,0.), (0.1,0.1,0.1), (5,10,4), api=3)
C._initVars(a, 'F', 1.); C._initVars(b, 'F', 2.); C._initVars(c, 'F', 3.)
C._initVars(a, 'centers:G', 1.); C._initVars(b, 'centers:G', 2.); C._initVars(c, 'centers:G', 3.)
a = T.join([a, b, c]); t = C.newPyTree(["Base", a])
#test.testT(t, 6)
