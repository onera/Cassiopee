# - splitBAR (array) -
import Transform as T
import Converter as C
import Generator as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (50,1,1))
a = C.convertArray2Tetra(a)
a = G.close(a)
b = T.splitBAR(a, 5)
test.testA(b, 1)
