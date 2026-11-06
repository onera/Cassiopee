# - ExtractMesh distribue -
import Converter.PyTree as C
import Distributor2.PyTree as Distributor2
import Converter.Mpi as Cmpi
import Transform.PyTree as T
import Post.Mpi as Pmpi
import Converter.Internal as Internal

rank = Cmpi.rank ; size = Cmpi.size

# Maillage d'extraction (distribue)
FILE = 'receiver.cgns'
e = Cmpi.convertFile2SkeletonTree(FILE)
(e, dic) = Distributor2.distribute(e, NProc=size, algorithm='fast', useCom=0)
e = Cmpi.readZones(e, FILE, rank=rank)

# Maillage source (distribue)
FILE = 'donor.cgns'
a = Cmpi.convertFile2SkeletonTree(FILE)
(a, dic) = Distributor2.distribute(a, NProc=size, algorithm='fast', useCom=0)
a = Cmpi.readZones(a, FILE, rank=rank)

# Extract mesh
e = Pmpi.extractMesh(a, e)

# Reconstruit l'arbre complet a l'ecriture
Cmpi.convertPyTree2File(e, 'out.cgns')
