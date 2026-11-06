# - center2Node distributed -
import Converter.PyTree as C
import Distributor2.PyTree as Distributor2
import Converter.Mpi as Cmpi
import Connector.PyTree as X
import Generator.PyTree as G

# Create test case
N = 11
t = C.newPyTree(['Base'])
pos = 0
for i in range(N):
    a = G.cart( (pos,0,0), (1,1,1), (10+i, 10, 10) )
    pos += 10 + i - 1
    t[2][1][2].append(a)
t = C.initVars(t, '{centers:Density} = {CoordinateX} + {CoordinateY}')
t = X.connectMatch(t)
if Cmpi.rank == 0: C.convertPyTree2File(t, 'in.cgns')
Cmpi.barrier()

# Reread in parallel
sk = Cmpi.convertFile2SkeletonTree('in.cgns')
(sk, dic) = Distributor2.distribute(sk, NProc=Cmpi.size, algorithm='gradient0',
                                    useCom='match')
a = Cmpi.readZones(sk, 'in.cgns', rank=Cmpi.rank)

# center2Node
a = Cmpi.center2Node(a, 'centers:Density')
# a is now a partial tree
a = C.rmVars(a, 'centers:Density')

# Rebuild full tree in file
Cmpi.convertPyTree2File(a, 'out.cgns')
