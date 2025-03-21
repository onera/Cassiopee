# - adapts a cells with respect to b points (array) -
import Intersector as XOR
import Converter as C
import Generator as G
import Converter.Internal as I
import numpy
import KCore.test as test

a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
aTH4 = C.convertArray2Tetra(a, split='withBarycenters')
a = C.convertArray2NGon(a); a = G.close(a)
aTH4 = C.convertArray2NGon(aTH4); aTH4 = G.close(aTH4)
#C.convertArrays2File([a], 'a.plt')

n = C.getNPts(a)
nodal_vals = numpy.empty((n,), dtype=I.E_NpyInt)
nodal_vals[:] = 2

## HEXA static adaptation
m = XOR.adaptCells(a, nodal_vals, sensor_type=2, smoothing_type=1)

m = XOR.closeCells(m[0])
#C.convertArrays2File([m], 'out.plt')
test.testA(m,1)

## TETRA static adaptation
n = C.getNPts(aTH4)
nodal_vals = numpy.empty((n,), dtype=I.E_NpyInt)
nodal_vals[:] = 2

m = XOR.adaptCells(aTH4, nodal_vals, sensor_type=2, smoothing_type=1)

m = XOR.closeCells(m[0])
#C.convertArrays2File([m], 'out1.plt')
test.testA(m,2)
