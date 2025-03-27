# - adapts a cells with respect to b points (array) -
import Intersector.PyTree as XOR
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as I
import Transform.PyTree as T
import Connector.PyTree as X
import numpy
import KCore.test as test

a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.convertArray2NGon(a); a = G.close(a)
#C.convertPyTree2File(a, 'a.cgns')
s = G.cartHexa((0.,0.,0.), (0.005,0.005,0.005), (5,5,5))

z = C.fillEmptyBCWith(a, 'wall', 'BCWall')
z = C.initVars(z, '{centers:Density} = {centers:CoordinateX} + {centers:CoordinateY}')

STYPE=3 # DIR
########################## create the hook
XOR._setZonesAndJoinsUId(z)
hmsh = XOR.createHMesh(z, subdiv_type=STYPE) # 0 : ISOTROPIC subdivision
########################################

#nodal specification
n = C.getNCells(z)
cell_vals = numpy.empty((n,), dtype=I.E_NpyInt)
cell_vals[:] = 1
# refine now with source mesh
z = XOR.adaptCells(z, cell_vals, subdiv_type=STYPE, sensor_type=3, itermax=-1, hmesh=hmsh)

z = XOR.conformizeHMesh(z, hmsh)     # each children faces replace its mother in any polyhedron

z = XOR.closeCells(z)            # close cells (adding point on lateral faces)

I._rmNodesByName(z, 'zid')
I._rmNodesByName(z, 'rid')

test.testT(z, 1)
#C.convertPyTree2File(z, "PT_t12.cgns")

########################## free the hook
XOR.deleteHMesh(hmsh)
#####################################
