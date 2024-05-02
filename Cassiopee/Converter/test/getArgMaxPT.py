# - getArgMax (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1.,1.,1.), (10,10,10))
argmax = C.getArgMax(a, 'CoordinateX'); print(argmax)
#>> [9.0]
