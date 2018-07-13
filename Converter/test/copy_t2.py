# - copy (array2) -
import Converter as C
import KCore.test as test

a = C.array('x,y,z',2,2,2,api=2)
b = C.copy(a)
test.testO(b, 1)
