# - silhouette (array) -
import Generator as G
import Converter as C
import Post as P
import Geom as D
import KCore.test as test


# structure
# 1D
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,1,1))
vector=[0.,1.,0.]
res = P.silhouette([a], vector)
test.testA(res,1)

vector=[1.,0.,0.]

# 2D
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,1,30))
a = C.initVars(a,'v',0.)
res = P.silhouette([a], vector)
test.testA(res,2)

# BAR
l1 = D.line((0.,0.,0.),(1.,0.,0.))
l2 = D.line((0.,0.,0.),(0.,1.,0.))
A = C.convertArray2Tetra([l1,l2]); A = C.initVars(A,'v',1.)
res =  P.silhouette(A,vector)
test.testA(res,3)

# TRI
A = C.convertArray2Tetra([a]); A = C.initVars(A,'v',1.)
res = P.silhouette(A, vector)
test.testA(res,4)

# QUAD
A = C.convertArray2Hexa([a]); A = C.initVars(A,'v',1.)
res = P.silhouette(A, vector)
test.testA(res,5)
