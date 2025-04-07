# - sharpEdges (array) -
import Converter as C
import Generator as G
import Post as P
import Transform as T
import Geom as D
import KCore.test as test

# structure
a1 = G.cart((0.,0.,0.),(1.5,1.,1.),(2,2,1))
a2 = T.rotate(a1,(0.,0.,0.),(0.,1.,0.),100.)
A = [a1,a2]; A = C.initVars(A,'v',0.)
res = P.sharpEdges(A,alphaRef=45.)
test.testA(res,1)

# QUAD
A = C.convertArray2Hexa([a1,a2]); A = C.initVars(A,'v',0.)
res = P.sharpEdges(A,alphaRef=45.)
test.testA(res,2)

# HEXA
A = C.convertArray2Tetra([a1,a2]); A = C.initVars(A,'v',0.)
res = P.sharpEdges(A, alphaRef=45.)
test.testA(res,3)

# 1D structure
l1 = D.line((0.,0.,0.),(1.,0.,0.))
l2 = D.line((0.,0.,0.),(0.,1.,0.))
res =  P.sharpEdges([l1,l2],alphaRef=45.)
test.testA(res,4)

# BAR
l1 = D.line((0.,0.,0.),(1.,0.,0.))
l2 = D.line((0.,0.,0.),(0.,1.,0.))
A = C.convertArray2Tetra([l1,l2]); A = C.initVars(A,'v',0.)
res =  P.sharpEdges(A,alphaRef=45.)
test.testA(res,5)

# NGON 1D
A = C.convertArray2NGon([l1,l2]); A = C.initVars(A,'v',0.)
res =  P.sharpEdges(A,alphaRef=45.)
test.testA(res,6)

# NGON 2D
A = C.convertArray2NGon([a1,a2]); A = C.initVars(A,'v',0.)
res = P.sharpEdges(A,alphaRef=45.)
test.testA(res,7)
