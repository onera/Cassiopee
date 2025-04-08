# - computeGraph (pyTree) -
# - match -
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Distributor2.PyTree as D2
import Generator.PyTree as G
import Connector.PyTree as X
import KCore.test as test

LOCAL = test.getLocal()

# Case
N = 11

# Cas test
t = C.newPyTree(['Base'])
off = 0
for i in range(N):
    a = G.cart( (off,0,0), (1,1,1), (10+i, 10, 10) )
    off += 9+i
    t[2][1][2].append(a)
t = X.connectMatch(t)
t, stats = D2.distribute(t, NProc=5, algorithm='fast', useCom=0)
if Cmpi.rank == 0: C.convertPyTree2File(t, LOCAL+'/in.cgns')
Cmpi.barrier()

# full
graph = Cmpi.computeGraph(t, type='match')
for i in graph:
    for k in graph[i]: graph[i][k].sort()
if Cmpi.rank == 0: test.testO(graph, 1)

# skel (doit etre identique)
t = Cmpi.convertFile2SkeletonTree(LOCAL+'/in.cgns')
graph = Cmpi.computeGraph(t, type='match')
for i in graph:
    for k in graph[i]: graph[i][k].sort()
if Cmpi.rank == 0: test.testO(graph, 2)

# load skel (doit etre identique)
t = Cmpi.convertFile2SkeletonTree(LOCAL+'/in.cgns')
t = Cmpi.readZones(t, LOCAL+'/in.cgns', rank=Cmpi.rank)
graph = Cmpi.computeGraph(t, type='match')
for i in graph:
    for k in graph[i]: graph[i][k].sort()
if Cmpi.rank == 0: test.testO(graph, 3)

# partial (doit etre identique)
t = Cmpi.convertFile2SkeletonTree(LOCAL+'/in.cgns')
procDict = D2.getProcDict(t)
t = Cmpi.readZones(t, LOCAL+'/in.cgns', rank=Cmpi.rank)
t = Cmpi.convert2PartialTree(t)
# procDict est requis pour les arbres partiels
# Le graph est correct uniquement quand on lance sur 5 procs
graph = Cmpi.computeGraph(t, type='match', procDict=procDict)
for i in graph:
    for k in graph[i]: graph[i][k].sort()
if Cmpi.rank == 0: test.testO(graph, 4)
