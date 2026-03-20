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

# standard array1 tests on closed curve
a = T.join([l1,l2,l3,l4])
a = T.consSmooth(a,sweeps)
test.testA([a], 1)

# standard array1 tests on open curve
b = T.join([l1,l2,l3])
b = T.consSmooth(b,sweeps)
C.convertArrays2File(b, 'out.plt')
test.testA([b], 2)

# standard bar1D tests on closed curve
c = C.convertArray2Tetra(a)
c = T.consSmooth(c,sweeps)
test.testA([c], 3)

# standard bar1D tests on open curve
d = C.convertArray2Tetra(b)
d = T.consSmooth(d,sweeps)
test.testA([d], 4)