# - getSmoothNormalMap (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test
a = G.cart((0.,0.,0.),(1.,1.,1.),(10,10,1))
b = G.cart((0.,0.,0.),(1.,1.,1.),(1,10,10))
b = T.rotate(b,(0.,0.,0.),(0.,1.,0.),45.)
#QUAD
c = C.convertArray2Hexa([a,b])
c = T.join(c); c = T.reorder(c,(1,))
c = T.rotate(c,(0.,0.,0.),(0.,1.,0.),15.)
c = G.getSmoothNormalMap(c,niter=4)
test.testT(c,1)
#QUAD
c = C.convertArray2Tetra([a,b])
c = T.join(c); c = T.reorder(c,(1,))
c = T.rotate(c,(0.,0.,0.),(0.,1.,0.),15.)
c = G.getSmoothNormalMap(c,niter=4)
test.testT(c,2)
