# - makeCartesian(pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test
a = G.cart((0.,0.,0.),(1.,1.,1.),(10,10,10))
a = T.reorder(a, (3,2,-1))
a = T.makeCartesianXYZ(a)
test.testT(a)
