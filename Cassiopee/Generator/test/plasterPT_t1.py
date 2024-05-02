# - plaster (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Transform.PyTree as T
import Post.PyTree as P
import Geom.PyTree as D
import KCore.test as test

a = D.sphere((0,0,0), 1, N=50)
a = T.subzone(a, (6,1,1), (25,100,1))
a = C.convertArray2Hexa(a); a = G.close(a)
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
# contours
c = P.exteriorFaces(a)
cs = T.splitConnexity(c)
cs = C.initVars(cs,'F',1.)

# plaster hole
p = G.plaster([cs[0]],[a])
test.testT(p)
