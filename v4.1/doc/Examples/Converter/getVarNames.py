# - getVarNames (array) -
import Converter as C
a = C.array('x,y,z,ro', 12, 9, 12)
print(C.getVarNames(a))
#>> ['x', 'y', 'z', 'ro']
