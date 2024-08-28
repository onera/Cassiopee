# - getIndexField (array) -
import Converter as C
import Generator as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,1,1))
n = C.getIndexField(a)
test.testA(n, 1)

a = G.cart((0,0,0), (1,1,1), (10,1,1), api=2)
n = C.getIndexField(a)
test.testA(n, 2)
