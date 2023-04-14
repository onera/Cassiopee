# - getSmoothNormalMap (array) -
import Converter as C
import Generator as G
import Transform as T
import KCore.test as test

a = G.cart((0.,0.,0.),(1.,1.,1.),(10,10,1))
b = G.cart((0.,0.,0.),(1.,1.,1.),(1,10,10))
b = T.rotate(b,(0.,0.,0.),(0.,1.,0.),45.)

# QUAD
c = C.convertArray2Hexa([a,b])
c = T.join(c); c = T.reorder(c,(1,))
c = T.rotate(c,(0.,0.,0.),(0.,1.,0.),15.)
s = G.getSmoothNormalMap(c, niter=4)
test.testA([s],1)

# TRI
c = C.convertArray2Tetra([a,b])
c = T.join(c); c = T.reorder(c,(1,))
c = T.rotate(c,(0.,0.,0.),(0.,1.,0.),15.)
s = G.getSmoothNormalMap(c, niter=4)
test.testA([s],2)

# STRUCT
c = T.join(a,b)
c = T.rotate(c,(0.,0.,0.),(0.,1.,0.),15.)
s = G.getSmoothNormalMap(c, niter=4)
test.testA([s],3)

# NGON
c = C.convertArray2NGon([a,b])
c = T.join(c); c = T.reorder(c,(1,))
c = T.rotate(c,(0.,0.,0.),(0.,1.,0.),15.)
s = G.getSmoothNormalMap(c, niter=4)
test.testA([s],4)
