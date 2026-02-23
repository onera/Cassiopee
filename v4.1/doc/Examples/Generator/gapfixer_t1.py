# - gapfixer (array) -
import Generator as G
import Converter as C
import Geom as D
import KCore.test as test

a = D.circle( (0,0,0), 1, N=100 )
a = C.convertArray2Tetra(a); a = G.close(a)
b = G.cart((-2.,-2.,0.), (0.1,0.1,1.), (50,50,1))
a1 = G.gapfixer(a, b)
test.testA([a1],1)

hp = D.point((0.5, 0.5, 0.))
a2 = G.gapfixer(a, b, hp, refine=0)
test.testA([a2],2)
