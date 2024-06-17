# - Calcul d'isosurfaces a la paroi en distribue -
import Converter.PyTree as C
import Distributor2.PyTree as Distributor2
import Converter.Mpi as Cmpi
import Transform.PyTree as T
import Post.PyTree as P

rank = Cmpi.rank ; size = Cmpi.size

# lecture du squelette
a = Cmpi.convertFile2SkeletonTree('cart.cgns')

# equilibrage
(a, dic) = Distributor2.distribute(a, NProc=size, algorithm='fast', useCom=0)

# load des zones locales dans le squelette
a = Cmpi.readZones(a, 'cart.cgns', rank=rank)

# Passage en arbre partiel
a = Cmpi.convert2PartialTree(a)

# Calcul de l'iso surface
isos = P.isoSurfMC(a, 'F', 0.3)
isos = Cmpi.setProc(isos, rank) # affecte les isos au proc local
t = C.newPyTree(['Base', isos])

# Reconstruit l'arbre complet a l'ecriture
Cmpi.convertPyTree2File(t, 'out.cgns')
