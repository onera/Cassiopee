# - extractMesh (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

# test arbre
# Maillage en noeuds
ni = 30; nj = 40; nk = 10
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m,'Density',1)
m = C.initVars(m,'MomentumX',1)
m = C.initVars(m,'centers:cellN',1)
m = C.addBC2Zone(m,'overlap','BCOverlap','imin')
m = C.fillEmptyBCWith(m, 'nref','BCFarfield')
t = C.newPyTree(['Base']); t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t[2][1][2].append(m)
a = G.cart((0.,0.,0.), (0.01, 0.01, 0.1), (40, 40, 1))
surf = P.extractMesh(t, a) 
test.testT(surf,1)

a = G.cart((0.,0.,0.), (0.01, 0.01, 0.1), (40, 40, 2))
a = C.addBC2Zone(a,'wall','BCWall','kmin')
a = C.fillEmptyBCWith(a,'ov','BCOverlap')
surf = P.extractMesh(t, a) 
test.testT(surf,2)
