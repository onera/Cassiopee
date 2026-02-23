# - blankIntersectingCells (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import Connector.PyTree as X
import KCore.test as test
import Transform.PyTree as T

body = D.sphere6((0,0,0),1.,10)
body = T.reorder(body,(-1,2,3))
dh = G.cart((0,0,0),(0.1,1,1),(10,1,1))
a = G.addNormalLayers(body,dh,niter=0)
C._initVars(a,'centers:cellN',1.)
C._initVars(a,'F',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'kmin')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'kmax' )
t = C.newPyTree(['Base',a])
t2 = X.blankIntersectingCells(t, depth=1)
test.testT(t2,1)
