# double-wall (para)
import Converter.PyTree as C
import Connector.PyTree as X
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import Generator.PyTree as G

import KCore.test as test

# full cylinder
a = G.cylinder((0,0,0), 0.5, 1., 0., 360, 1., (50,50,50)); a[0] = 'cyl1'
C._addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
bc = Internal.getNodeFromName(a, 'wall1')
C._tagWithFamily(bc, 'WALL1')
a = X.connectMatch(a)
C._fillEmptyBCWith(a, 'far', 'BCFarfield')
Cmpi._setProc(a, 0)

# small cylinder
b = G.cylinder((0,0,0.2), 0.5, 0.7, 20., 50, 0.5, (20,20,20)); b[0] = 'cyl2'
C._addBC2Zone(b, 'wall2', 'BCWall', 'jmin')
bc = Internal.getNodeFromName(b, 'wall2')
C._tagWithFamily(bc, 'WALL2')
C._fillEmptyBCWith(b, 'ov', 'BCOverlap')
Cmpi._setProc(b, 1)

t = C.newPyTree(['Base1', a, 'Base2', b])
if Cmpi.rank == 0: Internal._rmNodesByName1(t, 'Base2')
if Cmpi.rank == 1: Internal._rmNodesByName1(t, 'Base1')

C._initVars(t, 'centers:cellN', 1.)
X._applyBCOverlaps(t, depth=2, loc='centers')

Cmpi._addGhostCells(t, t, 2)

# without precond
tc = C.node2Center(t)
X._doubleWall(t, tc, 'WALL1', 'WALL2', ghostCells=True, check=False)
if Cmpi.rank == 1: test.testT(tc,1)

# with precond
tc = C.node2Center(t)
surfaceCenters1 = X.initDoubleWall(t, 'WALL1', check=False)
X._doubleWall(t, tc, 'WALL1', 'WALL2', ghostCells=True, check=False, surfaceCenters1=surfaceCenters1)
if Cmpi.rank == 1: test.testT(tc,2)
