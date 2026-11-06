# - projectOrthoSmooth (array) -
import Geom as D
import Generator as G
import Transform as T
import Converter as C
import KCore.test as test

# Structure
a = D.sphere((0,0,0), 1., 30)
b = G.cart((-0.5,-0.5,-1.5),(0.05,0.05,0.1), (20,20,1))
c = T.projectOrthoSmooth(b, [a], niter=2)
test.testA([c], 1)

# Element
a = D.sphere((0,0,0), 1., 30)
a = C.convertArray2Tetra(a)
b = G.cartHexa((-0.5,-0.5,-1.5),(0.05,0.05,0.1), (20,20,1))
c = T.projectOrthoSmooth(b, [a], niter=2)
test.testA([c], 2)

# NGON
a = D.sphere((0,0,0), 1., 30)
a = C.convertArray2NGon(a)
b = G.cartNGon((-0.5,-0.5,-1.5),(0.05,0.05,0.1), (20,20,1))
c = T.projectOrthoSmooth(b, [a], niter=2)
test.testA([c], 3)
