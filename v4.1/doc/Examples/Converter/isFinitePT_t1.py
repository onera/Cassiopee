# - isFinite (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._initVars(a, 'F', 1.)
test.testO(C.isFinite(a), 1)
test.testO(C.isFinite(a, var='CoordinateZ'), 2)
test.testO(C.isFinite(a, var='F'), 3)
