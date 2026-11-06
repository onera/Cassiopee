# - addProcNode (pyTree) -
import Generator.PyTree as G
import Distributor2.PyTree as D2
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = D2.addProcNode(a, 12)
test.testT(a, 1)
