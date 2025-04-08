# - getRegularityMap (array) -
import Generator as G
import Converter as C

a = G.cart( (0,0,0), (1,1,1), (50,50,1))
a = G.enforceX(a, 25, 0.1, 10, 10)
ac = C.node2Center(a)
reg = G.getRegularityMap(a)
reg = C.addVars([ac,  reg])
C.convertArrays2File([reg], "out.plt")
