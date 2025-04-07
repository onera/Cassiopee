# - getNCells (array2) -
import Converter as C
import KCore.test as test

a = C.array('x,y,z',2,2,2, api=2)
n = C.getNCells(a)
test.testO(n, 1)

a = C.array('x,y,z',12,4,'HEXA', api=2)
n = C.getNCells(a)
test.testO(n, 2)
