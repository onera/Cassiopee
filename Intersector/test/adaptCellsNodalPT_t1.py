# - adapts a cells with respect to b points (PyTree) -
import Intersector.PyTree as XOR
import Converter.PyTree as C
import Generator.PyTree as G
import numpy
import KCore.test as test

z = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
zTH4 = C.convertArray2Tetra(z, split='withBarycenters')
z = C.convertArray2NGon(z); z = G.close(z)
zTH4 = C.convertArray2NGon(zTH4); zTH4 = G.close(zTH4)
#C.convertPyTree2File([z], 'a.cgns')

n = C.getNPts(z)
#nodal_vals = numpy.zeros((n,), dtype=numpy.int32)
nodal_vals = numpy.empty((n,), dtype=numpy.int32)
nodal_vals[:] = 2

z = C.fillEmptyBCWith(z, 'wall', 'BCWall')

m = XOR.adaptCellsNodal(z, [nodal_vals])

m = XOR.closeOctalCells(m)
#C.convertPyTree2File(m, 'out.cgns')
test.testT(m,1)


# TETRA
zTH4 = C.fillEmptyBCWith(zTH4, 'wall', 'BCWall')
n = C.getNPts(zTH4)
#nodal_vals = numpy.zeros((n,), dtype=numpy.int32)
nodal_vals = numpy.empty((n,), dtype=numpy.int32)
nodal_vals[:] = 2

m = XOR.adaptCellsNodal(zTH4, [nodal_vals])

m = XOR.closeOctalCells(m)
#C.convertPyTree2File(m, 'out2.cgns')
test.testT(m,2)


