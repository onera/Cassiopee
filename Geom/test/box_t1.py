# - box (array) -
import Geom as D
import Converter as C
import KCore.test as test

a = D.box((0,0,0), (1,1,1))
test.testA(a, 1)
