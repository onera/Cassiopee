# - extractMesh (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

# test arbre
# Maillage en noeuds
ni = 30; nj = 40; nk = 10
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
C._initVars(m,'Density',1)
C._initVars(m,'MomentumX',1)
C._initVars(m,'centers:cellN',1)
C._addBC2Zone(m,'overlap','BCOverlap','imin')
C._fillEmptyBCWith(m, 'nref','BCFarfield')
t = C.newPyTree(['Base']); C._addState(t[2][1], 'Mach', 0.6)
t[2][1][2].append(m)
a = G.cart((0.,0.,0.), (0.01, 0.01, 0.1), (40, 40, 1))
P._extractMesh(t, a)
test.testT(a,1)

a = G.cart((0.,0.,0.), (0.01, 0.01, 0.1), (40, 40, 2))
C._addBC2Zone(a,'wall','BCWall','kmin')
C._fillEmptyBCWith(a,'ov','BCOverlap')
P._extractMesh(t, a)
test.testT(a,2)
