# - splitTBranches (pyTree)
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,1,50))
c1 = T.subzone(a,(1,1,1),(50,1,1))
c2 = T.subzone(a,(1,50,1),(50,50,1))
c3 = T.subzone(a,(1,1,1),(1,50,1))
c = [c1,c2,c3]; c = C.convertArray2Hexa(c)
c = T.join(c)
C._initVars(c,'F',1.); C._initVars(c,'centers:G',2.)
res = T.splitTBranches(c)
t = C.newPyTree(['Base', 1]); t[2][1][2] += res
test.testT(t)
