# - getStats (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1.,1.,1.), (11,1,1))
stats = C.getStats(a, 'x'); print(stats)
#>> {'min': 0.0, 'max': 10.0, 'mean': 5.0, 'stdv': 3.1622776601683795}