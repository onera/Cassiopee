# - getTriQualitylityMap (array) -
import Generator as G
import Converter as C
import Geom as D

a = D.sphere( (0,0,0), 1, N=10 )
a = C.convertArray2Tetra(a); a = G.close(a)
n = G.getTriQualityMap(a)
n = C.center2Node(n); n = C.addVars([a, n])
C.convertArrays2File([n], "out.plt")
