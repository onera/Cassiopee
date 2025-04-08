# - setIBCType (pyTree) -
import Converter.PyTree as C
import Geom.IBM as D_IBM
import Geom.PyTree as D
import KCore.test as test

a = D.circle((0,0,0), 1. , 0., 360.)
D_IBM._setIBCType(a, "Musker")
test.testT(a, 1)
