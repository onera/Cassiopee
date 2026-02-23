# - getDistantIndex -
import Converter as C
import Geom as D
import KCore.test as test

a = D.naca(12., 5001)
l = D.getLength(a)
l2 = D.getDistantIndex(a, 1, l/10.)
res = C.array('r', 1,1,1)
res[1][0,0] = l2

test.testA([res], 1)
