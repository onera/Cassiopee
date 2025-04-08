# - getRegularityMap (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Transform.PyTree as T
import KCore.test as test

# Test 3D structure
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,10))
a = C.addBC2Zone(a, 'wall1','BCWall','jmin')
a = C.addBC2Zone(a, 'wall2','BCWall','kmin')
a = C.addBC2Zone(a, 'match1','BCMatch','imin',a,'imax',[1,2,3])
a = C.addBC2Zone(a, 'match2','BCMatch','imax',a,'imin',[1,2,3])
a = C.fillEmptyBCWith(a, 'overlap', 'BCOverlap')
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
t = G.getRegularityMap(a)
test.testT(t, 1)

# Test 2D structure
msh = D.naca(12., 5001)
# Distribution
Ni = 300; Nj = 50
distrib = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))
a = G.hyper2D(msh, distrib, "C")
a = T.reorder(a, (-3,2,1))
a = C.addBC2Zone(a, 'wall1','BCWall','imin')
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
t = G.getRegularityMap(a)
test.testT(t, 2)

# Test 1D structure
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (1,10,1))
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
t = G.getRegularityMap(a)
test.testT(t, 3)
