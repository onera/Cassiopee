# - center2Node distribue -
import Converter.PyTree as C
import Distributor2.PyTree as Distributor2
import Converter.Mpi as Cmpi
import Connector.PyTree as X
import Generator.PyTree as G

# Case
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

# lecture du squelette
sk = Cmpi.convertFile2SkeletonTree('in.cgns')

# equilibrage
(sk, dic) = Distributor2.distribute(sk, NProc=Cmpi.size, algorithm='gradient0',
                                    useCom='match')

# load des zones locales dans le squelette
a = Cmpi.readZones(sk, 'in.cgns', proc=Cmpi.rank)

# center2Node
a = Cmpi.center2Node(a, 'centers:Density')
# a est maintenant un arbre partiel
a = C.rmVars(a, 'centers:Density')

# Reconstruit l'arbre complet a l'ecriture
Cmpi.convertPyTree2File(a, 'out.cgns')
