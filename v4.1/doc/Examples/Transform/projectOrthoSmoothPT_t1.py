# - projectOrthoSmooth (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

a = D.sphere((0,0,0), 1., 30)
b = G.cart((-0.5,-0.5,-1.5),(0.05,0.05,0.1), (20,20,1))
c = T.projectOrthoSmooth(b, [a], niter=2)
test.testT(c, 1)
