# Compute distance to walls prior to elsA computations

import Converter.PyTree as C
import Dist2Walls.PyTree as DTW
import Geom.PyTree as D
import Generator.PyTree as G
import Connector.PyTree as X

# Creation of the structured multiblock mesh
# kmin=1 define the wall BC
m = D.sphere6((0,0,0), 1., N=10)
di = G.cart((0.,0.,0.), (0.1,1.,1.), (11,1,1))
m = G.addNormalLayers(m, di)
m = C.addBC2Zone(m, 'wall', 'BCWall', 'kmin')
m = C.addBC2Zone(m, 'nref', 'BCFarfield', 'kmax')
m = X.connectMatch(m)
t = C.newPyTree(['Base']); t[2][1][2]=m
t = C.initVars(t, 'centers:cellN', 1.)
# 
# Get the wall surfaces 
#
walls = C.extractBCOfType(t, 'BCWall')
walls += C.extractBCOfType(t, 'BCWallViscous')
walls += C.extractBCOfType(t, 'BCWallInviscid')

# Convert the wall border to centers
bodies = C.node2Center(walls)
#
# Computes the distance field
#
t = DTW.distance2Walls(t, bodies, loc='centers', type='mininterf')
#
# Write the distance files for elsA
#
import numpy
import Converter.Internal as Internal
import Converter
nodes = Internal.getNodesFromName(t, 'TurbulentDistance')
c = 0
for n in nodes:
    ni = n[1].shape[0]; nj = n[1].shape[1]; nk = n[1].shape[2]
    a = numpy.reshape(n[1], (ni*nj*nk), order='Fortran')
    a = numpy.reshape(a, (1,ni*nj*nk))
    array = ['walldistance', a, ni, nj, nk]
    array = Converter.initVars(array, 'wallglobalindex', 1)
    Converter.convertArrays2File([array], 'dist'+str(c)+'.v3d',
                                 'bin_v3d')
    c += 1
