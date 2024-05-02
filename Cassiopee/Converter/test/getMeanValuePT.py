# - getMeanValue (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1.,1.,1.), (11,1,1))
meanval = C.getMeanValue(a, 'CoordinateX'); print(meanval)
#>> 5.0
