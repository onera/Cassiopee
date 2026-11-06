# - readNodesFromFilter (pyTree) -
import Converter.PyTree as C
import Converter.Filter as Filter
import Generator.PyTree as G

# Cree le fichier test
a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((12,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a,b])
C.convertPyTree2File(t, 'test.hdf')

# Relit les noeuds par un filtre
path = ['/Base/cart/GridCoordinates/CoordinateX',
        '/Base/cart/GridCoordinates/CoordinateY']

# start/stride/count (nbre d'entrees)/block
# Read 2x2x2 data from file, starting from 1,1,1
DataSpaceMMRY = [[0,0,0], [1,1,1], [2,2,2], [1,1,1]]
DataSpaceFILE = [[1,1,1], [1,1,1], [2,2,2], [1,1,1]]
DataSpaceGLOB = [[10,10,10]]

f = {}
f[path[0]] = DataSpaceMMRY+DataSpaceFILE+DataSpaceGLOB
f[path[1]] = DataSpaceMMRY+DataSpaceFILE+DataSpaceGLOB

# Lit seulement les chemins fournis, retourne un dictionnaire des chemins lus
a = Filter.readNodesFromFilter('test.hdf', f)
print(a[path[0]])
#>> [1. 2. 1. 2. 1. 2. 1. 2.]
