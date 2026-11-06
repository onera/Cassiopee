# - getMaxValue (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1.,1.,1.), (11,1,1))
maxval = C.getMaxValue(a, 'x'); print(maxval)
#>> 10.0
