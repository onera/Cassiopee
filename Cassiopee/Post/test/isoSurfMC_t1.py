# - isoSurfMC (array) -
import Post as P
import Converter as C
import Generator as G
import KCore.test as test

a = G.cartHexa((-20,-20,-20), (0.25,0.25,0.5), (100,100,50))
a = C.initVars(a, '{field}={x}*{x}+{y}*{y}+{z}')
iso = P.isoSurfMC(a, 'field', value=5.)
test.testA(iso, 1)
