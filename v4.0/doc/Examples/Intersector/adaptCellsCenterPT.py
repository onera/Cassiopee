# - adapts a cells with respect to b points (PyTree) -
import Intersector.PyTree as XOR
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import numpy

z = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
zTH4 = C.convertArray2Tetra(z, split='withBarycenters')
z = C.convertArray2NGon(z); z = G.close(z)
zTH4 = C.convertArray2NGon(zTH4); zTH4 = G.close(zTH4)
#C.convertPyTree2File([z], 'a.cgns')

n = C.getNCells(z)
nodal_vals = numpy.empty((n,), dtype=Internal.E_NpyInt)
nodal_vals[:] = 2

## HEXA static adaptation
z = C.fillEmptyBCWith(z, 'wall', 'BCWall')
z = C.initVars(z, '{centers:Density} = {centers:CoordinateX} + {centers:CoordinateY}')

XOR._setZonesAndJoinsUId(z)

m = XOR.adaptCells(z, [nodal_vals], sensor_type=3, smoothing_type=1)

m = XOR.closeCells(m)
C.convertPyTree2File(m, 'out.cgns')

## HEXA dynamic adaptation
hmsh = XOR.createHMesh(z, 0) # 0 : ISOTROPIC subdivision

m = XOR.adaptCells(z, [nodal_vals], sensor_type=3, hmesh=hmsh, smoothing_type=1)
m = XOR.closeCells(m)            # close cells (adding point on lateral faces)

C.convertPyTree2File(m, "out1.cgns")
XOR.deleteHMesh(hmsh);


## TETRA static adaptation
zTH4 = C.fillEmptyBCWith(zTH4, 'wall', 'BCWall')
zTH4 = C.initVars(zTH4, '{centers:Density} = {centers:CoordinateX} + {centers:CoordinateY}')

n = C.getNCells(zTH4)
nodal_vals = numpy.empty((n,), dtype=Internal.E_NpyInt)
nodal_vals[:] = 2

m = XOR.adaptCells(zTH4, [nodal_vals], sensor_type=3, smoothing_type=1)

m = XOR.closeCells(m)
C.convertPyTree2File(m, 'out2.cgns')
