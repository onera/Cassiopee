# - projectCloudSolution(pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import Post.PyTree as P
import Transform.PyTree as T
import Generator.PyTree as G
import KCore.test as test

a = D.sphere((0,0,0),1.,N=20)
a = C.convertArray2Tetra(a); a = G.close(a)
b = D.sphere6((0,0,0),1.,N=15)
b = C.convertArray2Tetra(b); b = T.join(b)
pts = C.convertArray2Node(b)
C._initVars(pts, '{F}={CoordinateX}*{CoordinateY}')
C._initVars(a, 'F', 0.)
P._projectCloudSolution(pts,a, loc='nodes')
test.testT(a)
P._projectCloudSolution(pts,a, loc='centers')
test.testT(a,2)
