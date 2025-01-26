# - getOrthogonalityMap (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as T

# Test 3D structure
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,10))
a = C.addBC2Zone(a, 'wall1','BCWall','jmin')
a = C.addBC2Zone(a, 'wall2','BCWall','kmin')
a = C.addBC2Zone(a, 'match1','BCMatch','imin',a,'imax',[1,2,3])
a = C.addBC2Zone(a, 'match2','BCMatch','imax',a,'imin',[1,2,3])
a = C.fillEmptyBCWith(a,'overlap','BCOverlap')
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
a = G.getOrthogonalityMap(a)
t = C.newPyTree(['Base']); t[2][1][2].append(a)
T.testT(t, 1)

# Test 2D structure
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a = C.addBC2Zone(a, 'wall1','BCWall','imin')
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
t = G.getOrthogonalityMap(a)
T.testT(t, 2)

# Test 2D structure
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (1,10,1))
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
t = G.getOrthogonalityMap(a)
T.testT(t, 3)
