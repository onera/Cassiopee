# - gapsmanager (array) -
import Geom as D
import Converter as C
import Generator as G
import Transform as T
import KCore.test as test

# Cas d'un maillage en centre
a = D.sphere6((0,0,0), 1, N=10)
a = C.node2Center(a)
a = C.convertArray2Tetra(a)
b = G.gapsmanager(a, mode=2)
test.testA(b, 1)

# Cas d'un maillage avec recouvrement
a = D.sphere6((0,0,0), 1, N=10)
a = C.convertArray2Tetra(a)
a = T.join(a)
a = G.close(a)
b = D.sphere6((0,0,0), 1, N=15)
b = b[0]
b = C.convertArray2Tetra(b)
b = T.rotate(b, (0,0,0), (0,0,1), 35.)
b = G.gapsmanager([a,b], mode=1)
test.testA(b, 2)
