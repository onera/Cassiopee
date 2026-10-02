# Usage: kpython -n 3 -t 4 connectMatchPT_m3.py
# - connectMatchNGon (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Post.PyTree as P
import Transform.PyTree as T
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Connector.Mpi as Xmpi
import KCore.test as test

# This test case is meant to replicate behavior
# that can be found in CODA simulations with overset.
# MPI processors are distributed over separate meshes
# with separate communicators. Based on the number of 
# meshes and procs, some meshes can be split while
# others are not. We need to make sure that the connectMatchNGon
# operation does not hang indefinitely in this condition.

N = 20

a1 = G.cartNGon((0.,0.,0.), (1./N,1./N,1), (N+1,N+1,1), api=3); a1[0] = 'zoneA'
#
a2 = T.translate(a1, (1,0,0)); a2[0] = 'zoneB'
#
a3 = G.cylinder((0.5,0.5,0), 0.1, 0.4, 0., 360., 0, (N,N,1))
a3 = C.convertArray2NGon(a3, api=3); a3[0] = 'zoneC'
p = P.exteriorFaces(a3)
p = T.splitConnexity(p)
C._addBC2Zone(a3, 'sym', 'BCSymmetryPlane', subzone=p[1])
C._addBC2Zone(a3, 'wall', 'BCWallViscous', subzone=p[0])

bc11 = G.cartNGon((0.,0.,0.), (1,1./N,1), (1,N+1,1), api=3); bc11[0] = 'bc11'
bc12 = G.cartNGon((0.,0.,0.), (1./N,1,1), (N+1,1,1), api=3); bc12[0] = 'bc12'
bc13 = G.cartNGon((0.,1.,0.), (1./N,1,1), (N+1,1,1), api=3); bc13[0] = 'bc13'

C._addBC2Zone(a1, 'farfield11', 'BCFarfield', subzone=bc11)
C._addBC2Zone(a1, 'farfield12', 'BCFarfield', subzone=bc12)
C._addBC2Zone(a1, 'farfield13', 'BCFarfield', subzone=bc13)

bc21 = T.translate(bc11, (2,0,0)); bc21[0] = 'bc21'
bc22 = T.translate(bc12, (1,0,0)); bc22[0] = 'bc22'
bc23 = T.translate(bc13, (1,0,0)); bc23[0] = 'bc23'

C._addBC2Zone(a2, 'farfield21', 'BCFarfield', subzone=bc21)
C._addBC2Zone(a2, 'farfield22', 'BCFarfield', subzone=bc22)
C._addBC2Zone(a2, 'farfield23', 'BCFarfield', subzone=bc23)

t = C.newPyTree(['Cart', [a1,a2], 'Cyl', [a3]])

for i, z in enumerate(Internal.getZones(t)):
    Cmpi._setProc(z, i)

# removes zones with proc != rank
Cmpi._convert2PartialTree(t, Cmpi.rank)
C._deleteEmptyBases(t)

# connectMatchNGon:
# The background grid is split in two.
# The curvilinear grid is not split and 
# all its exterior faces are already assigned a BC.
zones = Internal.getZones(t)
Xmpi._connectMatchNGon(zones[0])

if Cmpi.master: test.testT(t, 1)

