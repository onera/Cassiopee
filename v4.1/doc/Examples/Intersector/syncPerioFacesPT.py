# - concatenate (PyTree) -

import Converter.PyTree   as C
import Generator.PyTree   as G
import Connector.PyTree   as X
import Transform.PyTree   as T
import Intersector.PyTree as XOR
import Post.PyTree        as P
import Converter.Internal as I

## build a cartesian mesh
a    = G.cartHexa((0.,0.,0.), (1.,1.,1.), (10,10,10))
a    = C.convertArray2NGon(a)

XOR._reorient(a)

## triangulate top face #################
aF = P.exteriorFaces(a)
aF = T.splitSharpEdges(aF)
top = I.getZones(aF)[5]

BCs = [top]
BCNames = ['wall']
BCTypes = ['BCWallViscous']

C._recoverBCs(a, (BCs, BCNames, BCTypes))
a = XOR.triangulateBC(a, 'BCWallViscous')
C._deleteZoneBC__(a)
#########################################

a = XOR.syncMacthPeriodicFaces(a, translation=[0.,0.,9.])

a = X.connectMatchPeriodic(a, translation=[0.,0.,9.])

C.convertPyTree2File(a, 'out.cgns')
