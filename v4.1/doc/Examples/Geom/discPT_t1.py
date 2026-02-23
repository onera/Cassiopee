# - disc (pyTree) -
import Geom.PyTree as D
import KCore.test as test

a = D.disc((0,0,0), 1.)
test.testT(a, 1)
