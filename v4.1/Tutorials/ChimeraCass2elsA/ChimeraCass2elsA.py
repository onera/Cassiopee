# Chimera connectivity of 2 overlapping cylinders and a background grid
# Computation by Cassiopee and conversion to elsA files

import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.elsAProfile as CE
import numpy

# Chimera global parameter: 1 layer of interpolation cells
DEPTH=1

# Cylinder 1
a = G.cylinder((0,0,0),1.,3.,360,0,1,(200,30,2)); a[0] = 'cylindre1'
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'ov1', 'BCOverlap', 'jmax')

# Cylinder 2
b = G.cylinder((4,0,0),1.,3.,360,0,1,(200,30,2)); b[0] = 'cylindre2'
b = C.addBC2Zone(b, 'wall1', 'BCWall', 'jmin')
b = C.addBC2Zone(b, 'ov1', 'BCOverlap', 'jmax')

# Background grid
c = G.cart((-5.,-7.5,0), (15./200,15./200,1), (200,200,2))

# Creation of the pyTree
t = C.newPyTree(['Corps1', 'Corps2', 'Bgrd'])
t[2][1][2].append(a); t[2][2][2].append(b); t[2][3][2].append(c)

# BCs and joins
t = X.connectMatch(t, dim=2)
t = C.fillEmptyBCWith(t,'nref','BCFarfield', dim=2)

#==============================================================================
# 1. Set cellN to 2 to centers near overlap BCs
#==============================================================================
t = X.applyBCOverlaps(t, depth=DEPTH)

#==============================================================================
# 2. Compute blanking
#==============================================================================
# bodies description
bodyC1 = C.extractBCOfType(t[2][1][2], 'BCWall')
bodyC2 = C.extractBCOfType(t[2][2][2], 'BCWall')
bodies = [bodyC1,bodyC2]
# blanking matrix
BM = numpy.array([[0,1],[1,0],[1,1]])
# blanking
t = X.blankCells(t, bodies, BM, depth=DEPTH, dim=2)
# compute OversetHoles nodes
t = X.cellN2OversetHoles(t)
# set interpolated points around blanked points
t = X.setHoleInterpolatedPoints(t, depth=DEPTH)

#==============================================================================
# 3. Compute interpolation coefficients
#    and generate files readable by elsA with keys:
#    num.set('chm_conn_io','read')
#    num.set('chm_conn_fprefix','...')
#==============================================================================
# for interpolated cell centers
t = X.setInterpolations(t, loc='cell', prefixFile='interpolationFile',
                        sameBase=1)

# for interpolated interface (only for depth = 1)
if DEPTH == 1: t = X.setInterpolations(t, loc='face', prefixFile='interpolationFile', sameBase=1)

#==============================================================================
# 4. Optional: Convert tree to elsAxdt profile (for rereading by elsAxdt)
#==============================================================================
t = CE.convert2elsAxdt(t)

#==============================================================================
# 5. Save file
#==============================================================================
C.convertPyTree2File(t, 'out.cgns')
