# - enforceX (pyTree)-
import Generator.PyTree as G
import KCore.test as test
import Converter.PyTree as C

# 2D
Ni = 50; Nj = 50
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))
a = C.addBC2Zone(a, 'wall1','BCWall','jmin')
a = C.addBC2Zone(a, 'match1','BCMatch','imin',a,'imax',[1,2])
a = C.addBC2Zone(a, 'match2','BCMatch','imax',a,'imin',[1,2])
a = C.addBC2Zone(a, 'wall2','BCWall','jmax')
a = C.addVars(a,'Density'); a = C.addVars(a,'centers:cellN')
a2 = G.enforceX(a,0.3, 0.001, (13,25))
test.testT(a2,1)

# Exact number of added points
a2 = G.enforceX(a, 0.3, 0.001, 13,25)
test.testT(a2,2)

# 3D
Ni = 50; Nj = 50; Nk = 2
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,Nk))
a = C.addBC2Zone(a, 'overlap1','BCOverlap','imin')
a = C.addBC2Zone(a, 'overlap2','BCOverlap','imax')
a = C.addBC2Zone(a,'nref','BCFarfield',[1,1,5,10,1,2])
a = C.fillEmptyBCWith(a,'wall','BCWall',dim=3)
a = C.addVars(a,'Density'); a = C.addVars(a,'centers:cellN')
a2 = G.enforceX(a,0.3, 0.001, (13,25))
test.testT(a2,3)
