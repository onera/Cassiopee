# - map (array) -
import Geom as D
import Generator as G
import KCore.test as test
import Converter as C

l = D.line( (0,0,0), (1,1,0) )
Ni = 10
d = G.cart( (0,0,0), (1./(Ni-1),1.,1.), (Ni,1,1) )
m = G.map(l, d)
test.testA([m],1)

# Map on a curve
l = D.circle((0,0,0), 1. , 0., 360., 10)
Ni = 100
d = G.cart( (0,0,0), (1./(Ni-1),1.,1.), (Ni,1,1) )
l = C.convertArray2Tetra(l)
m = G.map(l, d)
test.testA([m],2)
