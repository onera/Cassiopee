# - projectCloudSolution (array) -
import Converter as C
import Geom as D
import Post as P
import Transform as T
import Generator as G
import KCore.test as test

a = D.sphere((0,0,0),1.,N=20)
a = C.convertArray2Tetra(a); a = G.close(a)
b = D.sphere6((0,0,0),1.,N=15)
b = C.convertArray2Tetra(b); b = T.join(b)
pts = C.convertArray2Node(b)
pts = C.initVars(pts,'{F}={x}*{y}')
a = C.initVars(a,'{F}=0.')
a1 = P.projectCloudSolution(pts, a)
test.testA([a1],1)
a = C.convertArray2Tetra(a); a = G.close(a)
a1 = P.projectCloudSolution(pts, a, loc='centers')
test.testA([a1],4)

a = D.naca(12., N=301)
pts = D.naca(12., N=351)
pts = C.convertArray2Node(pts)
pts = C.initVars(pts, '{F}={x}*{y}')
a = C.initVars(a, 'F', 0.)
a1 = P.projectCloudSolution(pts, a, dim=2)
test.testA([a1],2)
a = C.convertArray2Tetra(a); a = G.close(a)
a1 = P.projectCloudSolution(pts, a, dim=2, loc='centers')
test.testA([a1],3)
