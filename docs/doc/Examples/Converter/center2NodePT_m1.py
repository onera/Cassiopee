# - center2Node distribue -
import Converter.PyTree as C
import Distributor2.PyTree as Distributor2
import Converter.Mpi as Cmpi
import Connector.PyTree as X
import Generator.PyTree as G
import KCore.test as test

LOCAL = test.getLocal()

# Case
N = 11
t = C.newPyTree(['Base'])
pos = 0
for i in range(N):
    a = G.cart( (pos,0,0), (1,1,1), (10+i, 10, 10) )
    pos += 10 + i - 1
    t[2][1][2].append(a)
t = C.initVars(t, '{centers:Density} = {centers:CoordinateX} + {centers:CoordinateY}')
t = X.connectMatch(t)
if Cmpi.rank == 0: C.convertPyTree2File(t, LOCAL+'/in.cgns')
Cmpi.barrier()

# lecture du squelette
sk = Cmpi.convertFile2SkeletonTree(LOCAL+'/in.cgns')

# equilibrage
(sk, dic) = Distributor2.distribute(sk, NProc=Cmpi.size, algorithm='gradient0',
                                    useCom='match')

# load des zones locales dans le squelette
a = Cmpi.readZones(sk, LOCAL+'/in.cgns', rank=Cmpi.rank)

# center2Node
a = Cmpi.center2Node(a, 'centers:Density')
# a est maintenant un arbre partiel
a = C.rmVars(a, 'centers:Density')

if Cmpi.rank == 0: test.testT(a, 1)
Cmpi.barrier()
