# - getMassCenter (pyTree) -
import Geom.PyTree as D
import KCore.test as test

a = D.sphere((1,0,0), 1., N=30)
Xc = D.getMassCenter(a)
test.testO(Xc, 1)