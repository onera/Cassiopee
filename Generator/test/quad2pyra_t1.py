# - quad2Pyra (array) -
import Generator as G
import Converter as C
import Geom as D
import Transform as T
import KCore.test as test

a = D.sphere6((0,0,0), 1, N=10)
a = C.convertArray2Hexa(a)
a = T.join(a)
a = G.close(a)
b = G.quad2Pyra(a, hratio=1.)
test.testA([b],1)
