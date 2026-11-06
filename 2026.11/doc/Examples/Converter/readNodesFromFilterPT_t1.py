# - readNodesFromFilter (pyTree) -
import Converter.PyTree as C
import Converter.Filter as Filter
import Generator.PyTree as G
import KCore.test as test

LOCAL = test.getLocal()

# Cree le fichier test
a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cartHexa((12,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a,b])
C.convertPyTree2File(t, LOCAL+'/test.hdf')

# Relit les noeuds par un filtre
path = ['/Base/cart/GridCoordinates/CoordinateX',
        '/Base/cartHexa/GridCoordinates/CoordinateX']

f = {}
# start/stride/count (nbre d'entrees)/block
DataSpaceMMRY0 = [[0,0,0], [1,1,1], [2,2,2], [1,1,1]]
DataSpaceFILE0 = [[1,1,1], [1,1,1], [2,2,2], [1,1,1]]
DataSpaceGLOB0 = [[0]]
f[path[0]] = DataSpaceMMRY0+DataSpaceFILE0+DataSpaceGLOB0
DataSpaceMMRY1 = [[0], [1], [2], [1]]
DataSpaceFILE1 = [[0], [1], [2], [1]]
DataSpaceGLOB1 = [[0]]
f[path[1]] = DataSpaceMMRY1+DataSpaceFILE1+DataSpaceGLOB1

# Lit seulement les chemins fournis, retourne un dictionnaire des chemins lus
a = Filter.readNodesFromFilter(LOCAL+'/test.hdf', f)
test.testO(a)
#print a[path[0]].ravel('k')
#print a[path[1]].ravel('k')
