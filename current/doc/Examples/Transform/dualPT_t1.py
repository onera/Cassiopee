# - dual (pyTree)
import Converter.PyTree as C
import Transform.PyTree as T
import Geom.PyTree as D
import Generator.PyTree as G
import KCore.test as test

# NGON 2D ne contenant pas de faces externes
a = D.sphere((0,0,0), 1., 15)
a = C.convertArray2NGon(a); a = G.close(a)
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',1.);
res = T.dual(a)
t = C.newPyTree(['Base']); t[2][1][2]+=[res]
test.testT(t)

# NGON 2D contenant des faces externes
a = D.sphere((0,0,0), 1., 15)
a = T.subzone(a,(1,1,1),(15,15,1))
a = C.convertArray2NGon(a); a = G.close(a)
a = C.initVars(a,'F',1.)
res = T.dual(a)
t = C.newPyTree(['Base']); t[2][1][2] += [res]
test.testT(t,2)

# NGON 3D
a = G.cartHexa((0,0,0),(1,1,1),(11,11,11))
a = C.convertArray2NGon(a); a = G.close(a)
res = T.dual(a)
t = C.newPyTree(['Base']); t[2][1][2]+=[res]
test.testT(t,3)

a = G.cartTetra((0,0,0),(1,1,1),(11,11,11))
a = C.convertArray2NGon(a); a = G.close(a)
res = T.dual(a)
t = C.newPyTree(['Base']); t[2][1][2]+=[res]
test.testT(t,4)

a = D.sphere( (0,0,0), 1., 15); a = C.convertArray2Tetra(a); a = G.close(a)
distrib = G.cart((0,0,0),(0.1,0,0),(11,1,1))
a = G.addNormalLayers(a,distrib)
a = C.convertArray2NGon(a); a = G.close(a)
a = C.initVars(a,'F',1.)
res = T.dual(a)
t = C.newPyTree(['Base']); t[2][1][2]+=[res]
test.testT(t,5)

# NGON 1D - closed contour
a = D.circle((0,0,0), 1. , 0., 360., 10)
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',1.)
res = T.dual(a)
t = C.newPyTree(['Base']); t[2][1][2]+=[res]
test.testT(t,6)

# NGON 1D - open coutour
a = D.circle((0,0,0), 1. , 0., 180., 10)
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',1.)
res = T.dual(a)
t = C.newPyTree(['Base']); t[2][1][2]+=[res]
test.testT(t,7)
