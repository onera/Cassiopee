# - addChimera2Base (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cylinder((0,0,0), 1., 1.5, 0., 360., 1., (80,30,2))
a = C.addBC2Zone(a,'wall','BCWall','imin')
a = C.addBC2Zone(a,'overlap','BCOverlap','imax')
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)

t = C.newPyTree(['Base']); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t[2][1] = C.addChimera2Base(t[2][1], 'XRayTol', 1.e-6)
t[2][1] = C.addChimera2Base(t[2][1], 'XRayDelta', 0.1)
t[2][1] = C.addChimera2Base(t[2][1], 'DoubleWallTol', 100.)
t[2][1] = C.addChimera2Base(t[2][1], 'Priority', 1)
t[2][1] = C.addChimera2Base(t[2][1], '+', 'Base')
test.testT(t)
