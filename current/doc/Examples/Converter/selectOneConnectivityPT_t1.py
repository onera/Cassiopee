# - selectOneConnectivity (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Post.PyTree as P
import KCore.test as test

a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
p = P.exteriorFaces(a)
p = T.splitSharpEdges(p, 80.)
C._addBC2Zone(a, 'wall', 'BCWall', subzone=p[0])

# full
b = C.selectOneConnectivity(a, irange=[730,810])
test.testT(b, 1)
