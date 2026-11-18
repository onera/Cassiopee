# - getStats (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1.,1.,1.), (11,1,1))
stats = C.getStats(a, 'CoordinateX'); print(stats)
#>> {'min': 0.0, 'max': 10.0, 'mean': 5.0, 'stdv': 3.1622776601683795}