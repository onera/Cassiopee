# - getSmoothNormalMap (array) -
import Converter as C
import Generator as G
import Transform as T

a = G.cart((0.,0.,0.),(1.,1.,1.),(10,10,1))
b = G.cart((0.,0.,0.),(1.,1.,1.),(1,10,10))
b = T.rotate(b,(0.,0.,0.),(0.,1.,0.),45.)
c = C.convertArray2Hexa([a,b])
c = T.join(c); c = T.reorder(c,(1,))
c = T.rotate(c,(0.,0.,0.),(0.,1.,0.),15.)
s = G.getSmoothNormalMap(c, niter=4)
c = C.addVars([c,s])
C.convertArrays2File(c, "out.plt")
