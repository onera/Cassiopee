# - computeGraph (pyTree) -
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Distributor2.PyTree as Distributor2
import Generator.PyTree as G
import Connector.PyTree as X
import KCore.test as test

LOCAL = test.getLocal()

if Cmpi.rank == 0:
    a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
    b = G.cart( (9,0,0), (1,1,1), (10,10,10) )
    t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
    t = X.connectMatch(t)
    C.convertPyTree2File(t, LOCAL+'/in.cgns')
Cmpi.barrier()

t = Cmpi.convertFile2SkeletonTree(LOCAL+'/in.cgns')
(t, dic) = Distributor2.distribute(t, NProc=2, algorithm='fast')
t = Cmpi.readZones(t, LOCAL+'/in.cgns', rank=Cmpi.rank)

# Through BBox
tb = Cmpi.createBBoxTree(t)
graph = Cmpi.computeGraph(tb, type='bbox')
if Cmpi.rank == 0: test.testO(graph, 1)

# Through match
graph = Cmpi.computeGraph(t, type='match')
if Cmpi.rank == 0: test.testO(graph, 2)
