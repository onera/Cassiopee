# - getSharpestAngle (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Geom.PyTree as D

N = 10
d1 = G.cart((0.,0.,0.), (0.05,1,1),(N,1,4))
d2 = G.cart((0.,0.,0.), (1.,0.001,1),(1,10*N,4))
d2 = T.rotate(d2,(0.,0.,0.),(0.,0.,1.),30.)
s = T.join(d1,d2)
s = C.convertArray2Hexa(s)
s = T.reorder(s,(-1,))
s = D.getSharpestAngle(s)
C.convertPyTree2File(s, "out.cgns")
