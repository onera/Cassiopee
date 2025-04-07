# - breakConnectivity (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
b = G.cartTetra((9,0,0), (1,1,1), (10,10,10))
c = C.mergeConnectivity(a, b, boundary=0)

t = C.newPyTree(['Base',c])
t = C.breakConnectivity(t)
test.testT(t, 1)
