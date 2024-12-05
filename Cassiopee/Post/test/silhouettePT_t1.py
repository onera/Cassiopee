# - silhouette (array) -
import Generator.PyTree as G
import Converter.PyTree as C
import Post.PyTree as P
import Transform.PyTree as T
import Geom.PyTree as D
import KCore.test as test


# 1D structure
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,1,1))
a = C.initVars(a,'P',1.); a = C.initVars(a,'centers:Q',2.)
vector=[0.,1.,0.]
res = P.silhouette(a, vector)
test.testT(res,1)

vector=[1.,0.,0.]

# 2D
a1 = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,1,30))
a1 = C.addBC2Zone(a1,'wall','BCWall','jmin')
a1 = C.addBC2Zone(a1,'overlap','BCOverlap','jmax')
a2 = T.rotate(a1,(0.,0.,0.),(0.,1.,0.),120.); a2[0] = 'cart2'
A=[a1,a2]; A = C.initVars(A,'P',1.); A = C.initVars(A,'centers:Q',2.)
res = P.silhouette(A, vector)
test.testT(res,2)

# BAR
l1 = D.line((0.,0.,0.),(1.,0.,0.))
l2 = D.line((0.,0.,0.),(0.,1.,0.))
A = C.convertArray2Tetra([l1,l2]);
A = C.initVars(A,'P',1.); A = C.initVars(A,'centers:Q',2.)
res =  P.silhouette(A,vector)
test.testT(res,3)

# TRI
A = C.convertArray2Tetra([a1])
A = C.initVars(A,'P',1.); A = C.initVars(A,'centers:Q',2.)
res = P.silhouette(A, vector)
test.testT(res,4)

# QUAD
A = C.convertArray2Hexa([a1])
A = C.initVars(A,'P',1.); A = C.initVars(A,'centers:Q',2.)
res = P.silhouette(A, vector)
test.testT(res,5)
