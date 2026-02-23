# - refine (pyTree) -
import Geom.PyTree as D
import Transform.PyTree as T
import KCore.test as test

a = D.line((0,0,0), (1,0,0), N=10)
b = D.line((1,0,0), (2,1,0), N=30)
a = T.join([a,b])
a = D.refine(a, N=30)
test.testT(a, 1)
