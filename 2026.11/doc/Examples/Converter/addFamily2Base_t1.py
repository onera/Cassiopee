# - addFamily2Base (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cylinder((0,0,0), 1., 1.5, 0., 360., 1., (50,20,20))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'wall2', 'BCWall', 'imax')
t = C.newPyTree(['Base', a])
t[2][1] = C.addFamily2Base(t[2][1], 'flap', 'BCWall')
test.testT(t,1)
