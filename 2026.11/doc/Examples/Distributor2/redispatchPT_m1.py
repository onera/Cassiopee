# - redispatch (pyTree) -
import Converter.PyTree as C
import Distributor2.PyTree as D2
import Distributor2.Mpi as D2mpi
import Converter.Mpi as Cmpi
import Connector.PyTree as X
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

LOCAL = test.getLocal()

# Case
N = 11
t = C.newPyTree(['Base'])
pos = 0
for i in range(N):
    a = G.cart((pos,0,0), (1,1,1), (10+i, 10, 10))
    pos += 10 + i - 1
    t[2][1][2].append(a)
t = X.connectMatch(t)
if Cmpi.rank == 0: C.convertPyTree2File(t, LOCAL+'/in.cgns')
Cmpi.barrier()

# lecture du squelette
a = Cmpi.convertFile2SkeletonTree(LOCAL+'/in.cgns')

# equilibrage 1
(a, dic) = D2.distribute(a, NProc=Cmpi.size, algorithm='fast', useCom=0)

# load des zones locales dans le squelette
a = Cmpi.readZones(a, LOCAL+'/in.cgns', rank=Cmpi.rank)

# equilibrage 2 (a partir d'un squelette charge)
(a, dic) = D2.distribute(a, NProc=Cmpi.size, algorithm='gradient1',
                         useCom='match')
Cmpi._convert2PartialTree(a)
D2mpi._redispatch(a)
Internal._sortByName(a)

if Cmpi.rank == 0: test.testT(a, 1)

# force toutes les zones sur 0
zones = Internal.getZones(a)
for z in zones:
    node = Internal.getNodeFromName2(z, 'proc')
    Internal.setValue(node, 0)

D2mpi._redispatch(a)
Internal._sortByName(a)
if Cmpi.rank == 0: test.testT(a, 2)
