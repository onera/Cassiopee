# - getSharpestAngle (array) -
import Converter as C
import Generator as G
import Transform as T
import Geom as D
import KCore.test as test

N = 10
d1 = G.cart((0.,0.,0.), (0.05,1,1),(N,1,4))
d2 = G.cart((0.,0.,0.), (1.,0.001,1),(1,10*N,4))
d2 = T.rotate(d2,(0.,0.,0.),(0.,0.,1.),30.)
s0 = T.join(d1,d2)

# QUAD
s = C.convertArray2Hexa(s0)
s = T.reorder(s,(-1,))
r = D.getSharpestAngle(s)
s = C.addVars([s,r])
test.testA([s],1)

# TRI
s = C.convertArray2Tetra(s0)
s = T.reorder(s,(-1,))
r = D.getSharpestAngle(s)
s = C.addVars([s,r])
test.testA([s],2)

# BAR
d1 = G.cart((0.,0.,0.), (0.05,1,1),(N,1,1))
d2 = G.cart((0.,0.,0.), (0.05,0.001,1),(2*N,1,1))
d2 = T.rotate(d2,(0.,0.,0.),(0.,0.,1.),30.)
s0 = T.join(d1,d2)
s = C.convertArray2Hexa(s0)
r = D.getSharpestAngle(s)
s = C.addVars([s,r])
test.testA([s],3)
