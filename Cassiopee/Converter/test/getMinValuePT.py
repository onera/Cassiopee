# - getMinValue (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1.,1.,1.), (11,1,1))
C._initVars(a, '{centers:F}={centers:CoordinateX}')
minval = C.getMinValue(a, 'CoordinateX'); print(minval)
#>> 0.0
minval = C.getMinValue(a, 'centers:F'); print(minval)
#>> 0.5
minval = C.getMinValue(a, ['CoordinateX', 'CoordinateY']); print(minval)
#>> [0.0, 0.0]
minval = C.getMinValue(a, 'GridCoordinates'); print(minval)
#>> [0.0, 0.0, 0.0]
