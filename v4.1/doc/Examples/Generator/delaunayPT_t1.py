# - delaunay (pyTree)-
import Generator.PyTree as G
import KCore.test as test
import Converter.PyTree as C

# test structure 2d
ni = 11; nj = 11; nk = 1
hi = 1./(ni-1); hj = 1./(nj-1); hk = 1.
a = G.cart((0.,0.,0.), (hi,hj,hk), (20,1,20))
a = C.initVars(a,'centers:cellN',1); a = C.addVars(a,'Density')
b = G.delaunay(a)
test.testT(b,1)

# test non structure tri
ni = 11; nj = 11; nk = 1
hi = 1./(ni-1); hj = 1./(nj-1); hk = 1.
a = G.cartTetra((0.,0.,0.), (hi,hj,hk), (20,20,1))
a = C.initVars(a,'centers:cellN',1); a = C.addVars(a,'Density')
b = G.delaunay(a)
test.testT(b,2)

# test non structure quad
a = G.cartHexa((0.,0.,0.), (hi,hj,hk), (20,20,1))
b = G.delaunay(a)
test.testT(b,3)
