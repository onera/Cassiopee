# - setDfar (pyTree) -
import Converter.PyTree as C
import Geom.IBM as D_IBM
import Geom.PyTree as D
import KCore.test as test

a = D.circle((0,0,0), 1. , 0., 360.)
a = D_IBM.setDfar(a,10)
test.testT(a, 1)
