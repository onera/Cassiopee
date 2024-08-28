# - cartRx3 (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

a = G.cartRx3((0,0,0), (10,10,10), (0.1,0.1,0.1), (-10,-10,-10), (30,20,20), (1.1,1.1,1.1))

test.testT(a, 1)
