# - setValue (array) -
import Converter
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# test on array3
a = G.cart((0,0,0), (1,1,1), (10,10,10))
f = C.getFields('GridCoordinates', a, api=3)[0]
Converter.setValue(f, 12, (1,2,3))
ret = Converter.getValue(f, 12)
test.testO(ret, 1)
