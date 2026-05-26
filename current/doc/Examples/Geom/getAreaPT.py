# - getArea (pyTree) -
import Geom.PyTree as D
import math

a = D.sphere((0,0,0), 1., N=30)
area = D.getArea(a)
print(area, 4*math.pi*1*1)
