# - BB (pyTree) -
import Generator.PyTree as G
import Geom.PyTree as D
import KCore.test as test

s = D.circle((0,0,0), 1., N=100)
a = G.BB(s); a[0] = 'bbox'
test.testT(a, 1)
