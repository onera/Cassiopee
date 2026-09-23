# - center2Node (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.Mpi as Cmpi
import Connector.Mpi as Xmpi

import KCore.test as test

def F1(x,y,z): return x*y*z
def F2(x,y,z): return x*y

N = 5

# center2Node NGON + BCMatch
a = G.cartNGon((0.,0.,0.), (1./N,1./N,1./N), (N+1,N+1,N+1), api=3); a[0] = 'zoneA'
b = T.translate(a, (1,0,0)); b[0] = 'zoneB'
t = C.newPyTree(['Base', [a,b]])
C._initVars(t, 'centers:Star', 1.)
C._initVars(t, 'centers:Density', F1, ['centers:CoordinateX','centers:CoordinateY','centers:CoordinateZ'], isVectorized=True)
C._initVars(t, 'centers:Wars', 2.)
C._initVars(t, 'centers:Density2', F2, ['centers:CoordinateX','centers:CoordinateY','centers:CoordinateZ'], isVectorized=True)

for i, z in enumerate(Internal.getZones(t)):
    Cmpi._setProc(z, i)

# removes zones with proc != rank
Cmpi._convert2PartialTree(t, Cmpi.rank)

zones = Internal.getZones(t)
Xmpi._connectMatchNGon(zones[0])

t = Cmpi.center2Node(t, var=['centers:Density', 'centers:Density2'])
Internal._rmNodesByName(t, Internal.__FlowSolutionCenters__)

if Cmpi.master: test.testT(t, 1)

# center2Node NGON + BCMatch + BCDataset
a = G.cartNGon((0.,0.,0.), (1./N,1./N,1./N), (N+1,N+1,N+1), api=3); a[0] = 'zoneA'
b = T.translate(a, (1,0,0)); b[0] = 'zoneB'

bc1 = G.cartNGon((0.,0.,0.), (1,1./N,1./N), (1,N+1,N+1), api=3); bc1[0] = 'bc1'
bc2 = T.translate(bc1, (2,0,0)); bc2[0] = 'bc2'

C._addBC2Zone(a, 'inflow', 'BCInflow', subzone=bc1)
C._addBC2Zone(b, 'outflow', 'BCOutflow', subzone=bc2)

inflow = Internal.getNodeFromNameAndType(a, 'inflow', 'BC_t')
C._initBCDataSet(inflow, 'Density', 2)

outflow = Internal.getNodeFromNameAndType(b, 'outflow', 'BC_t')
C._initBCDataSet(outflow, 'Density', -2)

t = C.newPyTree(['Base', [a,b]])
C._initVars(t, 'centers:Density', F1, ['centers:CoordinateX','centers:CoordinateY','centers:CoordinateZ'], isVectorized=True)

for i, z in enumerate(Internal.getZones(t)):
    Cmpi._setProc(z, i)

# removes zones with proc != rank
Cmpi._convert2PartialTree(t, Cmpi.rank)

zones = Internal.getZones(t)
Xmpi._connectMatchNGon(zones[0])

t = Cmpi.center2Node(t, var=['centers:Density'])
Internal._rmNodesByName(t, Internal.__FlowSolutionCenters__)

if Cmpi.master: test.testT(t, 2)