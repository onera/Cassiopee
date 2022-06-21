# - line (pyTree) -
import Geom.IBM as IBM
import Geom.PyTree as D
import Converter.PyTree as C
import KCore.test as test

a = D.circle((0,0,0), 1. , 0., 360.)
a = IBM.setIBCType(a,"Musker")
test.testT(a, 1)
