# - sharpEdges (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import Transform.PyTree as T
import Geom.PyTree as D
import KCore.test as test

# liste de zones 2D structurees
a1 = G.cart((0.,0.,0.),(1.5,1.,1.),(2,2,1))
a1 = C.addBC2Zone(a1,'wall','BCWall','jmin')
a1 = C.addBC2Zone(a1,'overlap','BCOverlap','jmax')
a2 = T.rotate(a1,(0.,0.,0.),(0.,1.,0.),100.); a2[0] = 'cart2'
A = [a1,a2]; A = C.initVars(A,'P',1.); A = C.initVars(A,'centers:Q',2.)
res = P.sharpEdges(A,alphaRef=45.)
test.testT(res,1)

# liste de zones QUAD
A = C.convertArray2Hexa(A)
A = C.initVars(A,'P',1.); A = C.initVars(A,'centers:Q',2.)
res = P.sharpEdges(A,alphaRef=45.)
test.testT(res,2)

# liste de zones TRI
A = C.convertArray2Tetra(A)
A = C.initVars(A,'P',1.); A = C.initVars(A,'centers:Q',2.)
res = P.sharpEdges(A,alphaRef=45.)
test.testT(res,3)

# 1D structure
l1 = D.line((0.,0.,0.),(1.,0.,0.))
l2 = D.line((0.,0.,0.),(0.,1.,0.))
A = [l1,l2]; A = C.initVars(A,'P',1.); A = C.initVars(A,'centers:Q',2.)
res =  P.sharpEdges(A,alphaRef=45.)
test.testT(res,4)

# BAR
l1 = D.line((0.,0.,0.),(1.,0.,0.))
l2 = D.line((0.,0.,0.),(0.,1.,0.))
A = C.convertArray2Tetra([l1,l2])
A = C.initVars(A,'P',1.); A = C.initVars(A,'centers:Q',2.)
res =  P.sharpEdges(A,alphaRef=45.)
test.testT(res,5)
