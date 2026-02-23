# - getVarNames (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
a = G.cart((0,0,0),(1,1,1),(10,10,10))
C._addVars(a, ['Density', 'centers:cellN'])
# one zone
print(C.getVarNames(a, loc='nodes'))
#>> [['CoordinateX', 'CoordinateY', 'CoordinateZ', 'Density']]
print(C.getVarNames(a, loc='centers'))
#>> [['CoordinateX', 'CoordinateY', 'CoordinateZ', 'centers:cellN']]
print(C.getVarNames(a, excludeXYZ=True, loc='both'))
#>> [['Density', 'centers:cellN']]
