# - quad2Pyra (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Transform.PyTree as T
import KCore.test as test

a = D.sphere6((0,0,0), 1, N=10)
a = C.convertArray2Hexa(a)
a = T.join(a)
a = G.close(a)
b = G.quad2Pyra(a, hratio=1.)
test.testT(b,1)
