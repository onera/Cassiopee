# Chimera assembly of 2 overlapping cylinders and a background Cartesian grid 

import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import numpy

# Cylinder 1
a = G.cylinder((0,0,0),1,10,360,0,1,(200,30,2)); a[0] = 'cylindre1'
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
# Cylinder 2
b =  G.cylinder((4,0,0),1,10,360,0,1,(200,30,2)); b[0] = 'cylindre2'
b = C.addBC2Zone(b, 'wall1', 'BCWall', 'jmin')
# Background grid
c = G.cart((-15.,-15.,0.),(0.3,0.3,1.),(111,111,2))
c = C.fillEmptyBCWith(c,'nref','BCFarfield',dim=2)
# Creation of the pyTree
t = C.newPyTree(['Corps1','Corps2','Bgrd'])
t[2][1][2].append(a); t[2][2][2].append(b); t[2][3][2].append(c)
# BCs
t = X.connectMatch(t, dim=2)
t = C.fillEmptyBCWith(t, 'overlap', 'BCOverlap', dim=2)
C.convertPyTree2File(t, 'basic0.cgns')

#==============================================================================
# 1. Set cellN to 2 to centers near overlap BCs
#==============================================================================
t = X.applyBCOverlaps(t, depth=2)
#==============================================================================
# 2. Blanking
#==============================================================================
# bodies description
bodyC1 = C.extractBCOfType(t[2][1][2], 'BCWall')
bodyC2 = C.extractBCOfType(t[2][2][2], 'BCWall')
bodies = [bodyC1,bodyC2]
# blanking matrix
BM = numpy.array([[0,1],[1,0],[1,1]])
# blanking
t = X.blankCells(t, bodies, BM, depth=2, dim=2)
# set interpolated cells around interpolated points
t = X.setHoleInterpolatedPoints(t, depth=2)
#==============================================================================
# 3. Overlap optimisation with a high priority to cylinders
#==============================================================================
t = X.optimizeOverlap(t, priorities=['Corps1',0,'Corps2',0,'Bgrd',1])
t = X.maximizeBlankedCells(t, depth=2)
#==============================================================================
# 4. Convert indices of cellN = 0 into OversetHoles nodes
#==============================================================================
t = X.cellN2OversetHoles(t)
#==============================================================================
# 5. Computes Chimera connectivity
#==============================================================================
t = X.setInterpolations(t, loc='cell')
#==============================================================================
# 6. Save file
#==============================================================================
C.convertPyTree2File(t, 'comp0.cgns')
