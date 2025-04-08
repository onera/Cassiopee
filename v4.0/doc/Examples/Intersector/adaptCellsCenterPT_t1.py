# - adapts a cells with respect to b points (PyTree) -
import Intersector.PyTree as XOR
import Converter.PyTree as C
import Converter.Internal as I
import Generator.PyTree as G
import numpy
import KCore.test as test

z = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
zTH4 = C.convertArray2Tetra(z, split='withBarycenters')
z = C.convertArray2NGon(z); z = G.close(z)
zTH4 = C.convertArray2NGon(zTH4); zTH4 = G.close(zTH4)
#C.convertPyTree2File([z], 'a.cgns')

n = C.getNCells(z)
nodal_vals = numpy.empty((n,), dtype=I.E_NpyInt)
nodal_vals[:] = 2

## HEXA static adaptation
z = C.fillEmptyBCWith(z, 'wall', 'BCWall')
z = C.initVars(z, '{centers:Density} = {centers:CoordinateX} + {centers:CoordinateY}')

XOR._setZonesAndJoinsUId(z)

m = XOR.adaptCells(z, [nodal_vals], sensor_type=3, smoothing_type=1, subdiv_type=0, hmesh=None)

m = XOR.closeCells(m)

I._rmNodesByName(m, 'zid')
I._rmNodesByName(m, 'rid')

#C.convertPyTree2File(m, 'nodalPT_t1_1.cgns')
test.testT(m,1)

## HEXA dynamic adaptation
hmsh = XOR.createHMesh(z, 0) # 0 : ISOTROPIC subdivision

m = XOR.adaptCells(z, [nodal_vals], sensor_type=3, smoothing_type=1, subdiv_type=0, hmesh=hmsh)
m = XOR.conformizeHMesh(m, hmsh)
m = XOR.closeCells(m)            # close cells (adding point on lateral faces)

I._rmNodesByName(m, 'zid')
I._rmNodesByName(m, 'rid')

XOR.deleteHMesh(hmsh);
#C.convertPyTree2File(m, 'nodalPT_t1_2.cgns')
test.testT(m,2)



## TETRA static adaptation
zTH4 = C.fillEmptyBCWith(zTH4, 'wall', 'BCWall')
zTH4 = C.initVars(zTH4, '{centers:Density} = {centers:CoordinateX} + {centers:CoordinateY}')

n = C.getNCells(zTH4)
nodal_vals = numpy.empty((n,), dtype=I.E_NpyInt)
nodal_vals[:] = 2

XOR._setZonesAndJoinsUId(zTH4)

m = XOR.adaptCells(zTH4, [nodal_vals], sensor_type=3, smoothing_type=1, subdiv_type=0)

m = XOR.closeCells(m)

I._rmNodesByName(m, 'zid')
I._rmNodesByName(m, 'rid')

#C.convertPyTree2File(m, 'nodalPT_t1_3.cgns')
test.testT(m,3)
