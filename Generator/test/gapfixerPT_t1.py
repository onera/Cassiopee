# - gapfixer (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

a = D.circle((0,0,0),1,N=100)
a = C.convertArray2Tetra(a); a = G.close(a)
b = G.cart((-2.,-2.,0.), (0.1,0.1,1.), (50,50,1))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a1 = G.gapfixer(a, b)
test.testT(a1, 1)

hp = D.point((0.5, 0.5, 0.))
a2 = G.gapfixer(a, b, hp, refine=0)
test.testT(a2, 2)
