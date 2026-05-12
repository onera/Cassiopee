# - getArea (pyTree) -
import Geom.PyTree as D
import KCore.test as test

a = D.sphere((0,0,0), 1., N=30)
area = D.getArea(a)
test.testO(area, 1)
