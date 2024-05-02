# - normL0 (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (11,11,11))
C._initVars(a, 'centers:F', 1.)
print('normL0 =', C.normL0(a, 'centers:F'))
#>> normL0 = 1.0
