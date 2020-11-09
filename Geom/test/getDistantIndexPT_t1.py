# - getDistantIndex (PyTree)-
import Geom.PyTree as D
import KCore.test as test

a = D.naca(12., 5001)
l = D.getLength(a)
l2 = D.getDistantIndex(a, 1, l/10.)
test.testO(l2, 1)
