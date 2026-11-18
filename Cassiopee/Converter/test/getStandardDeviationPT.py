# - getStandardDeviation (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1.,1.,1.), (11,1,1))
stdval = C.getStandardDeviation(a, 'CoordinateX'); print(stdval)
#>> 3.1622776601683795
C._initVars(a, '{centers:F}={centers:CoordinateX}')
stdval = C.getStandardDeviation(a, 'centers:F'); print(stdval)
#>> 2.8722813232690143