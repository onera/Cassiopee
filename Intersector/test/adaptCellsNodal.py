# - adapts a cells with respect to b points (array) -
import Intersector as XOR
import Converter as C
import Generator as G
import numpy

a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.convertArray2NGon(a); a = G.close(a)
#C.convertArrays2File([a], 'a.plt')

n = C.getNPts(a)
nodal_vals = numpy.empty((n,), dtype=numpy.int32)
nodal_vals[:] = 2

m = XOR.adaptCellsNodal(a, nodal_vals)

m = XOR.closeCells(m[0])
C.convertArrays2File([m], 'out.plt')


