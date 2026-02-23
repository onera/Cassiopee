# - splitBAR (PyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Transform.PyTree as T
import KCore.test as test

a = G.cart( (0,0,0), (1,1,1), (50,1,1) )
a = C.convertArray2Tetra(a)
a = G.close(a)
B = T.splitBAR(a, 5)
test.testT(B, 1)
