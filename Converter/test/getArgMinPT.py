# - getArgMin (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1.,1.,1.), (10,10,10))
argmin = C.getArgMin(a, 'CoordinateX'); print(argmin)
#>> [0.0]
