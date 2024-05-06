# - PR2VL -

import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import Converter.converter
import Connector.PyTree as X

# Extrait les champs vertex d'une condition aux limites
a = G.cart((0,0,0), (1,1,1), (10,10,1))
C._addBC2Zone(a, 'wall','BCWall', 'jmin')
n = Internal.getNodeFromName(a, 'PointRange')[1]
inds = Converter.converter.PR2VL(n, 10, 10, 1)
array2 = C.getFields('coords', a, api=2)[0]
sub = Converter.converter.extractFields(array2, inds)

a = G.cart((0,0,0), (1,1,1), (10,10,1)); a[0] = 'cart1'
b = G.cart((9,0,0), (1,1,1), (10,10,1)); b[0] = 'cart2'
C._initVars(a, '{ux}=0.')
C._initVars(a, '{uy}=0')
C._initVars(a, '{uz}=0')
C._initVars(b, '{ux}=1.')
C._initVars(b, '{uy}=1')
C._initVars(b, '{uz}=1')

t = C.newPyTree(['Base',a,b])
t = X.connectMatch(t, dim=2)

# Recupere un raccord
r = Internal.getNodeFromType(a, 'GridConnectivity1to1_t')
n = Internal.getNodeFromName(r, 'PointRange')[1]
nd = Internal.getNodeFromName(r, 'PointRangeDonor')[1]
trf = Internal.getNodeFromName(r, 'Transform')[1]

# Recupere les indices des raccords
inds, donors = Converter.converter.PR2VL(n, 10, 10, 1, nd, trf, 10, 10, 1)
print(inds)
print(donors)

array1 = C.getFields('nodes', a, vars=['ux'], api=2)[0]
array2 = C.getFields('nodes', b, vars=['ux'], api=2)[0]

# Recupere les champs correspondant aux indices
sub1 = Converter.converter.extractFields(array1, inds)
sub2 = Converter.converter.extractFields(array2, donors)
#print(sub2)

# Reimporte
Converter.converter._setPartialFieldsAverage(array1, [inds], [sub2])
print(array1)

C.convertPyTree2File(t, 'out.cgns')
