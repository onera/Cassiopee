# - map (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

l = D.line((0,0,0),(1,1,0))
l = C.initVars(l,'Density',1.); l = C.initVars(l,'centers:cellN',1.)
Ni = 11; dist = G.cart( (0,0,0), (1./(Ni-1),1.,1.), (Ni,1,1) )
l2 = G.map(l, dist)
test.testT(l2)

# Map on a curve
l = D.circle((0,0,0), 1. , 0., 360., 10)
Ni = 100
d = G.cart( (0,0,0), (1./(Ni-1),1.,1.), (Ni,1,1) )
l = C.convertArray2Tetra(l)
l = C.initVars(l,'Density',1.); l = C.initVars(l,'centers:cellN',1.)
m = G.map(l, d)
test.testT(m,2)

# 2D
a = G.cart((0,0,0),(1,1,1),(11,11,1))
a = C.addBC2Zone(a, 'wall', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'overlap', 'BCOverlap', [11,11,1,5,1,1])
a = C.fillEmptyBCWith(a, 'nref', 'BCFarfield', dim=2)
Ni = 15; dist = G.cart((0,0,0), (1./(Ni-1),1.,1.), (Ni,1,1))
a2 = G.map(a, dist, dir=2)
test.testT(a2, 3)
