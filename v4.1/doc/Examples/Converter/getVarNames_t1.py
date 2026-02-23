# - getVarNames (array) -
import Converter as C
import KCore.test as test

# Array
a = C.array('x,y,z,ro', 12, 9, 12)
b = C.array('x,F', 13, 2, 1)
val = C.getVarNames(a)
test.testO(val, 1)

# Arrays
val = C.getVarNames([a,b])
test.testO(val, 2)
