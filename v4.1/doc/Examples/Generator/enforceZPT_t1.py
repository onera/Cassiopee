# - enforceZ (pyTree)-
import Generator.PyTree as G
import KCore.test as test
import Converter.PyTree as C

Ni = 50; Nj = 2; Nk = 50
a = G.cart((0,0,0), (1./(Ni-1), 1., 0.5/(Nk-1)), (Ni,Nj,Nk))
a = C.addBC2Zone(a, 'wall1','BCWall','jmin')
a = C.addBC2Zone(a, 'overlap1','BCOverlap','imin')
a = C.addBC2Zone(a, 'overlap2','BCOverlap','imax')
a = C.addBC2Zone(a, 'wall2','BCWall','jmax')
a = C.addVars(a,'Density'); a = C.addVars(a,'centers:cellN')

# Exact number of added points
b = G.enforceZ(a, 0.3, 0.001, 10,15)
test.testT(b,1)
