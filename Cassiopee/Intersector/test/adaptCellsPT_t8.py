import Converter.PyTree as C
import Converter.Internal as I
import Generator.PyTree as G
#import Post.PyTree as P
import os, sys
import Intersector.PyTree as XOR
import numpy
import KCore.test as test

a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.convertArray2NGon(a); a = G.close(a)
#C.convertPyTree2File(a, 'a.cgns')
s = G.cartHexa((0.,0.,0.), (0.005,0.005,0.005), (5,5,5))

z = C.fillEmptyBCWith(a, 'wall', 'BCWall')
z = C.initVars(z, '{centers:Density} = {centers:CoordinateX} + {centers:CoordinateY}')


########################## create the hook
XOR._setZonesAndJoinsUId(z)
hmsh = XOR.createHMesh(z, 0) # 0 : ISOTROPIC subdivision
########################################

for i in range(3): # simple loop to demonstrate the feature sequencing adaptCells and adaptCellsNodal

    #nodal specification
    n = C.getNPts(z)
    nodal_vals = numpy.empty((n,), dtype=I.E_NpyInt)
    nodal_vals[:] = 1
    #one nodal_vals and one hmesh per zone
    z = XOR.adaptCellsNodal(z, [nodal_vals], hmesh=hmsh)
    # refine now with source mesh
    z = XOR.adaptCells(z, s, itermax=-1, hmesh=hmsh)

z = XOR.conformizeHMesh(z, hmsh)     # each children faces replace its mother in any polyhedron

z = XOR.closeCells(z)            # close cells (adding point on lateral faces)

I._rmNodesByName(z, 'zid')
I._rmNodesByName(z, 'rid')

test.testT(z, 1)
#C.convertPyTree2File(z, "PT_t8.cgns")

########################## free the hook
XOR.deleteHMesh(hmsh);
#####################################
