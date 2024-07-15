# - Calcul de distance a la paroi en distribue -
import Converter.PyTree as C
import Distributor2.PyTree as Distributor2
import Converter.Mpi as Cmpi
import Transform.PyTree as T
import Dist2Walls.PyTree as Dist2Walls

rank = Cmpi.rank ; size = Cmpi.size

# lecture des parois
bodies = C.convertFile2PyTree('walls.cgns')

# lecture du squelette
a = Cmpi.convertFile2SkeletonTree('in.cgns')

# equilibrage
(a, dic) = Distributor2.distribute(a, NProc=size, algorithm='fast', useCom=0)

# load des zones locales dans le squelette
a = Cmpi.readZones(a, 'in.cgns', rank=rank)

# Passage en arbre partiel
a = Cmpi.convert2PartialTree(a)

# Calcul des distances a la paroi
a = Dist2Walls.distance2Walls(a, bodies, type='mininterf')

# Reconstruit l'arbre complet a l'ecriture
Cmpi.convertPyTree2File(a, 'out.cgns')
