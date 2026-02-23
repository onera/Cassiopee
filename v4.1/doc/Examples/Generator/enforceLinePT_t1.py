# - enforceLine (pyTree) -
import Generator.PyTree as G
import KCore.test as test
import Geom.PyTree as D
import Converter.PyTree as C

# Distribution
Ni = 50; Nj = 50
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))
b = D.line((0.,0.2,0.), (1.,0.2,0.), 20)
a = C.addBC2Zone(a, 'wall1','BCWall','jmin')
a = C.addBC2Zone(a, 'match1','BCMatch','imin',a,'imax',[1,2])
a = C.addBC2Zone(a, 'match2','BCMatch','imax',a,'imin',[1,2])
a = C.addBC2Zone(a, 'wall2','BCWall','jmax')
a = C.addVars(a,'Density'); a = C.initVars(a,'centers:cellN',1)
c = G.enforceLine(a, b, 0.01, (5,3))
test.testT(c,1)
