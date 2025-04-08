# - Blanking en distribue -
import Converter.PyTree as C
import Distributor2.PyTree as Distributor2
import Converter.Mpi as Cmpi
import Transform.PyTree as T
import Connector.PyTree as X
import Converter.Internal as Internal
import numpy

rank = Cmpi.rank ; size = Cmpi.size

# lecture des corps servant a masquer
bodies = C.convertFile2PyTree('walls.cgns')
bodies = Internal.getNodesFromType(bodies, 'Zone_t')

# lecture du squelette
a = Cmpi.convertFile2SkeletonTree('in.cgns')

# equilibrage
(a, dic) = Distributor2.distribute(a, NProc=size, algorithm='fast', useCom=0)

# load des zones locales dans le squelette
a = Cmpi.readZones(a, 'in.cgns', rank=rank)

# Passage en arbre partiel
a = Cmpi.convert2PartialTree(a)

# Blanking local
BM = numpy.array([[1]])
a = X.blankCells(a, [bodies], BM, blankingType='node_in', dim=3)

# Reconstruit l'arbre complet a l'ecriture
Cmpi.convertPyTree2File(a, 'out.cgns')
