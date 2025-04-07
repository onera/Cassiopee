# - center2Node (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1,1,1), (30,40,1))
a = C.initVars(a, 'ro', 1.)
an = C.center2Node(a)
C.convertArrays2File(an, "out.plt")
