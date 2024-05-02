# - getMeanValue (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1.,1.,1.), (11,1,1))
meanval = C.getMeanValue(a, 'x'); print(meanval)
#>> 5.0
