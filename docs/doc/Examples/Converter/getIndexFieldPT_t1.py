# - getIndexField (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,1,1))
C._getIndexField(a)
test.testT(a, 1)
