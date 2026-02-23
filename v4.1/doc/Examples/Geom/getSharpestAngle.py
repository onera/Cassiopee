# - getSharpestAngle (array) -
import Converter as C
import Generator as G
import Transform as T
import Geom as D

N = 10
d1 = G.cart((0.,0.,0.), (0.05,1,1),(N,1,4))
d2 = G.cart((0.,0.,0.), (1.,0.001,1),(1,10*N,4))
d2 = T.rotate(d2,(0.,0.,0.),(0.,0.,1.),30.)
s = T.join(d1,d2)
s = C.convertArray2Hexa(s)
s = T.reorder(s,(-1,))
r = D.getSharpestAngle(s)
s = C.addVars([s,r])
C.convertArrays2File(s, "out.plt")
