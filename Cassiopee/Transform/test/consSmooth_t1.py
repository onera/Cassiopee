# - consSmooth (array) -
import Transform as T
import Geom as D
import KCore.test as test
import Converter as C

sweeps = 5
l1 = D.line((0.,0.,0.), (0.,1.,0.), N=5)
l2 = D.line((0.,1.,0.), (1.,1.,0.), N=5)
l3 = D.line((1.,1.,0.), (1.,0.,0.), N=5)
l4 = D.line((1.,0.,0.), (0.,0.,0.), N=5)

# struct 1D closed curve
a = T.join([l1,l2,l3,l4])
a = T.consSmooth(a, sweeps)
test.testA([a], 1)

# struct 1D open curve
b = T.join([l1,l2,l3])
b = T.consSmooth(b, sweeps)
test.testA([b], 2)

# bar1D closed curve
c = C.convertArray2Tetra(a)
c = T.consSmooth(c, sweeps)
test.testA([c], 3)

# bar1D open curve
d = C.convertArray2Tetra(b)
d = T.consSmooth(d, sweeps)
test.testA([d], 4)