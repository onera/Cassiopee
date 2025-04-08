# - extractMesh (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

# Maillage en noeuds
ni = 30; nj = 40; nk = 10
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m,'Density',1)
m = C.initVars(m,'MomentumX',1)
m = C.initVars(m,'centers:cellN',1)

a = G.cart((0.,0.,0.), (0.01, 0.01, 0.1), (40, 40, 1))
for i in [2]:
    a2 = P.extractMesh(m, a, order=i)
    test.testT(a2, i)
