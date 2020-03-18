# - adapts a cells with respect to b points (array) -
import Intersector as XOR
import Converter as C
import Generator as G
import numpy
import KCore.test as test

a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
aTH4 = C.convertArray2Tetra(a, split='withBarycenters')
a = C.convertArray2NGon(a); a = G.close(a)
aTH4 = C.convertArray2NGon(aTH4); aTH4 = G.close(aTH4)
#C.convertArrays2File([a], 'a.plt')

n = C.getNPts(a)
nodal_vals = numpy.empty((n,), dtype=numpy.int32)
nodal_vals[:] = 2

## HEXA static adaptation
m = XOR.adaptCellsNodal(a, nodal_vals)

m = XOR.closeCells(m[0])
test.testA(m,1)

## TETRA static adaptation
n = C.getNPts(aTH4)
nodal_vals = numpy.empty((n,), dtype=numpy.int32)
nodal_vals[:] = 2

m = XOR.adaptCellsNodal(aTH4, nodal_vals)

m = XOR.closeCells(m[0])
test.testA(m,2)
