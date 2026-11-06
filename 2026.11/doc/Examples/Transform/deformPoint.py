# - deformPoint (array) -
import Generator as G
import Transform as T
import Converter as C

a1 = G.cart((0,0,0), (1,1,1), (10,10,1))
a2 = T.deformPoint(a1, (0,0,0), (0.1,0.1,0.1), 0.5, 2.)
C.convertArrays2File([a2], "out.plt")
