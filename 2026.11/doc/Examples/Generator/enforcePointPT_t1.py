# - enforcePoint (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

Ni = 20; Nj = 20
a = G.cart((0,0,0), (1./(Ni-1),5./(Nj-1),1), (Ni,Nj,1))
a = C.addBC2Zone(a, 'wall1','BCWall','jmin')
a = C.addBC2Zone(a, 'match1','BCMatch','imin',a,'imax',[1,2])
a = C.addBC2Zone(a, 'match2','BCMatch','imax',a,'imin',[1,2])
a = C.addBC2Zone(a, 'wall2','BCWall','jmax')
a = C.addVars(a,'Density'); a = C.addVars(a,'centers:cellN')
b = G.enforcePoint(a, 0.5)
test.testT(b,1)
