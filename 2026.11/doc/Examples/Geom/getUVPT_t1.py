# - getUV (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.PyTree as C
import Post.PyTree as P
import KCore.test as test

a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
a = P.exteriorFaces(a)
a = C.initVars(a, 'u', 0.)
a = C.initVars(a, 'v', 0.)
ret = D.getUV(a, 2., 1920)

test.testT(ret[0], 0)
