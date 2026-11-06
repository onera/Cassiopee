# - isFinite (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._initVars(a, 'F', 1.)
print(C.isFinite(a))
#>> True
print(C.isFinite(a, var='CoordinateZ'))
#>> True
print(C.isFinite(a, var='F'))
#>> True
