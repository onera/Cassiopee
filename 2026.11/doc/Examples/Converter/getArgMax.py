# - getArgMax (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1.,1.,1.), (10,10,10))
argmax = C.getArgMax(a, 'x'); print(argmax)
#>> [9.0, 0.0, 0.0]
