# - setFluidInside (pyTree) -
import Converter.PyTree as C
import Geom.IBM as D_IBM
import Geom.PyTree as D
import KCore.test as test

#Prepare Geometry
a = D.circle((0,0,0), 1. , 0., 360.)

D_IBM._setFluidInside(a)
test.testT(a)
