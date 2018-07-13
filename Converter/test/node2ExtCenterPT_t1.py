# - node2ExtCenter (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

def F(x,y): return 2*x+y

# maillage 2D
ni = 30; nj = 40; nk = 1
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
a = C.initVars(a,'F',F,['CoordinateX','CoordinateY'])
ac = C.node2ExtCenter(a)
t = C.newPyTree(['Base',2]); t[2][1][2].append(ac)
test.testT(t,1)

# maillage 3D
ni = 30; nj = 40; nk = 30
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),10./(nk-1)), (ni,nj,nk))
a = C.initVars(a,'F',F,['CoordinateX','CoordinateY'])
ac = C.node2ExtCenter(a)
t = C.newPyTree(['Base']); t[2][1][2].append(ac)
test.testT(t,2)
