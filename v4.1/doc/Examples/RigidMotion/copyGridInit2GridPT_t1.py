# - copyGridInit2Grid (pyTree) -
import Generator.PyTree as G
import RigidMotion.PyTree as R
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (5,5,5))
R._copyGrid2GridInit(a, mode=1)
R._copyGridInit2Grid(a)
test.testO(a, 1)
