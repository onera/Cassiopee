# - isNamePresent (PyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (50,50,50))
C._initVars(a, 'F', 1.)
C._initVars(a, 'centers:G', 0.)

b = G.cart((0,0,0), (1,1,1), (50,50,50))
C._initVars(b, 'F', 2.)
C._initVars(b, 'centers:H', 3.)

t = C.newPyTree(['Base',a,b])

print(C.getVarNames([a, b]))
#>> [['CoordinateX', 'CoordinateY', 'CoordinateZ', 'F', 'centers:G'], ['CoordinateX', 'CoordinateY', 'CoordinateZ', 'F', 'centers:H']]

print(C.isNamePresent(a, 'F'))
#>> 1
print(C.isNamePresent(a, 'centers:F'))
#>> -1
print(C.isNamePresent([a, b], 'F'))
#>> 1
print(C.isNamePresent([a, b], 'centers:G'))
#>> 0
