# - getArea (pyTree) -
import Geom.PyTree as D
import numpy

a = D.sphere((0,0,0), 1., N=30)
area = D.getArea(a)
print(area, 4*numpy.pi*1*1)
