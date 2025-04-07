# - fittingPlaster (array) -
import Generator as G
import Converter as C
import Geom as D
import KCore.test as test

a = D.circle( (0,0,0), 1, N=50 )
a = C.convertArray2Tetra(a)
a = G.close(a)
b = G.fittingPlaster(a, bumpFactor=0.5)
test.testA([b])
