# - addPointInDistribution (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# Distribution
Ni = 50; Nj = 50
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,2))
a = C.addBC2Zone(a, 'match1','BCMatch','imin',a,'imax',[1,2])
a = C.addBC2Zone(a, 'match2','BCMatch','imax',a,'imin',[1,2])
a = C.addBC2Zone(a, 'wall2','BCWall','jmax')
a = C.addVars(a,'Density'); a = C.initVars(a,'centers:cellN',1)
b = G.addPointInDistribution( a, Ni )
t = C.newPyTree(['Base',2]); t[2][1][2].append(b)
test.testT(t,1)
