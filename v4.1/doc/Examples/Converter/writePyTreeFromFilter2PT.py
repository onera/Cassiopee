# - writePyTreeFromFilter (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Filter as Filter
import Converter.Internal as Internal

# Cas test
a = G.cart((0,0,0),(1,1,1),(10,10,10))
t = C.newPyTree(['Base', a])

# Filter
DataSpaceMMRY = [[0,0,0], [1,1,1], [10,10,10], [1,1,1]]
DataSpaceFILE = [[0,0,0], [1,1,1], [10,10,10], [1,1,1]]
DataSpaceGLOB = [[10,10,10]]
fr = {}
fr['/Base/cart/GridCoordinates/CoordinateX'] = DataSpaceMMRY+DataSpaceFILE+DataSpaceGLOB

# Formation squelette pour le filtre
ts = Internal.copyRef(t)
for path in fr:
    node = Internal.getNodeFromPath(ts, path)
    node[1] = None

# Ecriture squelette sur le proc 0
C.convertPyTree2File(ts, 'out.hdf')

# Dimensionnement des tableaux du filtre sur le proc 0 (skelData=None)
Filter.writePyTreeFromFilter(t, 'out.hdf', fr, skelData=None)

# Ecriture de la data partielle effective (skelData=[])
Filter.writePyTreeFromFilter(t, 'out.hdf', fr, skelData=[])

# Si HDF parallele, nous ne sommes pas oblige de dimensionner.
