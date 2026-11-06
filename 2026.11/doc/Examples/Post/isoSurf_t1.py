# - isoSurf (array) -
import Post as P
import Converter as C
import Generator as G
import Geom as D
import KCore.test as test

# Volume tetra
a = G.cartTetra((-20,-20,-20), (0.25,0.25,0.5), (100,100,50))
a = C.initVars(a, '{field}={x}*{x}+{y}*{y}+{z}')
iso = P.isoSurf(a, 'field', value=5.)
test.testA(iso, 1)

# Surface
b = D.sphere( (0,0,0), 1.)
b = C.initVars(b, '{field}={x}*{x}+{y}*{y}+{z}')
iso = P.isoSurf(b, 'field', value=0.5)
test.testA(iso, 2)

# Melting pot
c = G.cart((-10,-10,-10), (1,1,1), (10,10,10))
c = C.initVars(c, '{field}={x}*{x}+{y}*{y}+{z}')
all = [a,b,c]
iso = P.isoSurf(all, 'field', value=0.5)
test.testA(iso, 3)
