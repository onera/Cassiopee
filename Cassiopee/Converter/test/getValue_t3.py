# - getValue (array) -
import Converter as C
import Generator.PyTree as GP
import Converter.PyTree as CP
import KCore.test as test

# test on array3
a = GP.cart((0,0,0), (1,1,1), (10,10,10))
f = CP.getFields('GridCoordinates', a, api=3)[0]
ret = C.getValue(f, 12)
#C.setValue(f, 12, ret)
test.testO(ret, 1)
