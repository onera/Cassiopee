# - center2Node (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import Post.PyTree as P
import Transform.PyTree as T

import KCore.test as test

N = 5

# center2Node NGON + BCDataset (multi fields)
a = G.cartNGon((0.,0.,0.), (1./N,1./N,1./N), (N+1,N+1,N+1), api=3); a[0] = 'zoneA'

p = P.exteriorFaces(a)
C._addBC2Zone(a, 'farfield', 'BCFarfield', subzone=p)

bc = Internal.getNodeFromNameAndType(a, 'farfield', 'BC_t')
C._initBCDataSet(bc, 'Density', 12)
C._initBCDataSet(bc, 'VelocityX', 100)

t = C.newPyTree(['Base', a])
C._initVars(t, 'centers:Density', 0)
C._initVars(t, 'centers:VelocityX', 1)
C._initVars(t, 'centers:VelocityY', 2)

t = C.center2Node(t, var=['centers:Density', 'centers:VelocityX', 'centers:VelocityY'])
Internal._rmNodesByName(t, Internal.__FlowSolutionCenters__)

test.testT(t, 1)

# center2Node NGON + BCDataset (multi BCs)
a = G.cartNGon((0.,0.,0.), (1./N,1./N,1./N), (N+1,N+1,N+1), api=3); a[0] = 'zoneA'

p = P.exteriorFaces(a)
p = T.splitSharpEdges(p, 75.)
for i, pLoc in enumerate(p):
    C._addBC2Zone(a, 'farfield%d'%i, 'BCFarfield', subzone=pLoc)
    bc = Internal.getNodeFromNameAndType(a, 'farfield%d'%i, 'BC_t')
    C._initBCDataSet(bc, 'Density', i)

t = C.newPyTree(['Base', a])
C._initVars(t, 'centers:Density', 0)

t = C.center2Node(t, var='centers:Density')
Internal._rmNodesByName(t, Internal.__FlowSolutionCenters__)

test.testT(t, 2)
