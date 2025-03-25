# - normL2 (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (11,11,11))
a = C.initVars(a, "F", 1.)
print('normL2 =', C.normL2(a, "F"))
#>> normL2 = 1.0
