# - getSharpestAngle (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Geom.PyTree as D
import KCore.test as test

N = 10
d1 = G.cart((0.,0.,0.), (0.05,1,1),(N,1,4))
d2 = G.cart((0.,0.,0.), (1.,0.001,1),(1,10*N,4))
d2 = T.rotate(d2,(0.,0.,0.),(0.,0.,1.),30.)
s0 = T.join(d1,d2)
# QUAD
s = C.convertArray2Hexa(s0)
s = T.reorder(s,(-1,))
s = D.getSharpestAngle(s)
test.testT(s,1)
#TRI
s = C.convertArray2Tetra(s0)
s = T.reorder(s,(-1,))
s = D.getSharpestAngle(s)
test.testT(s,2)
