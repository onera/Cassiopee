# - writePyTreeFromFilter (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import Converter.Filter as Filter

# Ecrit la zone en entier
t = C.newPyTree(['Base'])
a = G.cart((0,0,0), (1,1,1), (10,10,10))
t[2][1][2] += [a]
C.convertPyTree2File(t, 'out.hdf')

# Prend une subzone et la remplace dans le fichier
t = C.newPyTree(['Base'])
b = T.subzone(a, (2,2,2), (5,5,5)); b[0] = 'cart'
C._initVars(b, 'CoordinateX', 1.)
t[2][1][2] += [b]

DataSpaceMMRY = [[0,0,0], [1,1,1], [4,4,4], [1,1,1]]
DataSpaceFILE = [[2,2,2], [1,1,1], [4,4,4], [1,1,1]]
DataSpaceGLOB = [[10,10,10]]

f = {}
f['/Base/cart/GridCoordinates/CoordinateX'] = DataSpaceMMRY+DataSpaceFILE+DataSpaceGLOB

# skelData != None since node already exists
Filter.writePyTreeFromFilter(t, 'out.hdf', f, skelData=[])
