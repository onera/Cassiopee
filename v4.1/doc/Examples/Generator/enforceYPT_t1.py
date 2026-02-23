# - enforceY (pyTree)-
import Generator.PyTree as G
import KCore.test as test
import Converter.PyTree as C

# 2D
Ni = 50; Nj = 50
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))
a = C.addBC2Zone(a, 'wall1','BCWall','jmin')
a = C.addBC2Zone(a, 'overlap1','BCOverlap','imin')
a = C.addBC2Zone(a, 'overlap2','BCOverlap','imax')
a = C.addBC2Zone(a, 'wall2','BCWall','jmax')
a = C.addVars(a,'Density'); a = C.addVars(a,'centers:cellN')
a2 = G.enforceY(a,0.3, 0.001, (10,15))
test.testT(a2,1)

# Exact number of added points
a2 = G.enforceY(a, 0.3, 0.001, 10,15)
test.testT(a2,2)
