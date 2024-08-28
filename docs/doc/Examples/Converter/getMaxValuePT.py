# - getMaxValue (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1.,1.,1.), (11,2,2))
C._initVars(a, '{centers:F}={centers:CoordinateX}')
maxval = C.getMaxValue(a, 'CoordinateX'); print(maxval)
#>> 10.0
maxval = C.getMaxValue(a, 'centers:F'); print(maxval)
#>> 9.5
maxval = C.getMaxValue(a, ['CoordinateX', 'CoordinateY']); print(maxval)
#>> [10.0, 1.0]
maxval = C.getMaxValue(a, 'GridCoordinates'); print(maxval)
#>> [10.0, 1.0, 1.0]
