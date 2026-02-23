# - stitchedHat (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test
import Converter.PyTree as C

c = D.circle( (0,0,0), 1., 360., 0., 100)
c = T.contract(c, (0,0,0), (0,1,0), (0,0,1), 0.1)
c = C.initVars(c,'centers:cellN',1.); c = C.initVars(c,'Density', 2.)
c = G.stitchedHat(c, (0,0,0), 1.e-4)
t = C.newPyTree(['Base',2]); t[2][1][2].append(c)
test.testT(t,1)
