# - cartRx (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# cartRx avec depth=1 et cellN
a = G.cartRx((0,0,0), (1,1,1), (10,10,10), (5,5,5), depth=1, addCellN=True)
test.testT(a, 1)