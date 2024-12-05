# - getTriQualityMap (array) -
import Geom as D
import Generator as G
import Converter as C
import KCore.test as test

a = D.sphere((0,0,0), 1, 50)
a = C.convertArray2Tetra(a)
a = G.close(a)
n = G.getTriQualityMap(a)
test.testA([n],1)
