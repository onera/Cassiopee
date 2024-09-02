# - getEdgeRatio(array) -
import Generator as G
import Converter as C

a = G.cart((0.,0.,0.), (0.1,0.1,0.1), (11,11,11))
a = G.enforcePlusX(a,1e-6,(5,50))
r = G.getEdgeRatio(a)
r = C.center2Node(r); r = C.addVars([a,r])
C.convertArrays2File([r], "out.plt")
