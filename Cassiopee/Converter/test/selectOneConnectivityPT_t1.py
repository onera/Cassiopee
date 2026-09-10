# - selectOneConnectivity (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import Transform.PyTree as T
import Post.PyTree as P
import KCore.test as test

# -- BE
# full, by irange
a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
a = C.initVars(a, 'F', 0.)
a = C.initVars(a, 'centers:G', 1.)
p = P.exteriorFaces(a)
p = T.splitSharpEdges(p, 80.)
C._addBC2Zone(a, 'wall', 'BCWall', subzone=p[0])
t = Internal.copyTree(a)
b = C.selectOneConnectivity(t, irange=[730,810])
test.testT(b, 1)

# full, by irange, keep volumic
t = Internal.copyTree(a)
b = C.selectOneConnectivity(t, irange=[1,729])
test.testT(b, 2)

# full, by name
a = G.cartHexa((0,0,0), (1,1,1), (5,5,5))
sz1 = G.cartHexa((0,0,0), (1,1,0), (5,5,1)); sz1[0] = "QUAD1"
sz2 = G.cartHexa((0,0,4), (1,1,1), (5,5,1)); sz2[0] = "QUAD2"
C._addBC2Zone(a, 'wall1', 'BCWall', subzone=sz1)
C._addBC2Zone(a, 'wall2', 'BCWall', subzone=sz2)
a = C.initVars(a, 'F', 0.)
a = C.initVars(a, 'centers:G', 1.)
t = Internal.copyTree(a)
b = C.selectOneConnectivity(t, name='QUAD2')
test.testT(b, 3)

# full, by number
t = Internal.copyTree(a)
b = C.selectOneConnectivity(t, number=1)
test.testT(b, 4)

# slice by irange, volumic
t = Internal.copyTree(a)
b = C.selectOneConnectivity(t, irange=[2,7])
test.testT(b, 5)

# slice by irange, surfacic
t = Internal.copyTree(a)
b = C.selectOneConnectivity(t, irange=[81,93])
test.testT(b, 6)

# wrong inputs
t = Internal.copyTree(a)
b = C.selectOneConnectivity(t, name='QUAD3')
test.testT(b, 7)

t = Internal.copyTree(a)
b = C.selectOneConnectivity(t, number=10)
test.testT(b, 8)

t = Internal.copyTree(a)
b = C.selectOneConnectivity(t, irange=[200,230])
test.testT(b, 9)

t = Internal.copyTree(a)
b = C.selectOneConnectivity(t, irange=[1,64+1])
test.testT(b, 10)

t = Internal.copyTree(a)
b = C.selectOneConnectivity(t, irange=[65,93])
test.testT(b, 11)

# -- NGON
