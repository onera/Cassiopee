# - getMinValue (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1.,1.,1.), (11,1,1))
minval = C.getMinValue(a, 'x'); print(minval)
#>> 0.0
