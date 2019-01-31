# - addGhostCells (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Transform.PyTree as T
import Connector.PyTree as X
import KCore.test as test

Nj = 5
a = G.cart((1,1,1), (1.,1.,1.), (4,Nj,6)); a[0]='cart1'
b = G.cart((4,1,1), (1.,1.,1.), (3,Nj,3)); b[0]='cart2'
c = G.cart((4,1,3), (1.,1.,1.), (3,Nj,6)); c[0]='cart3'
d = G.cart((1,1,6), (1.,1.,1.), (4,Nj,3)); d[0]='cart4'
a = T.reorder(a,(-2,1,3))
b = T.reorder(b,(1,2,3))
c = T.reorder(c,(3,1,2))

# init cellN
a = C.initVars(a, 'centers:cellN', 1.)
b = C.initVars(b, 'centers:cellN', 0.)
c = C.initVars(c, 'centers:cellN', 1.)
d = C.initVars(d, 'centers:cellN', 0.)

t = C.newPyTree(['Base',a,b,c,d])
t = C.initVars(t, '{F}=3*{CoordinateX}+2*{CoordinateY}')
t = C.initVars(t, '{centers:G}={centers:CoordinateZ}')
t = X.connectMatch(t,dim=3)
t = C.fillEmptyBCWith(t,'wall','BCWall')
t = Internal.addGhostCells(t,t,2,adaptBCs=0,fillCorner=1)
test.testT(t,1)

t = C.newPyTree(['Base',a,b,c,d])
t = C.initVars(t, '{F}=3*{CoordinateX}+2*{CoordinateY}')
t = C.initVars(t, '{centers:G}={centers:CoordinateZ}')
t = X.connectMatch(t,dim=3)
t = C.fillEmptyBCWith(t,"wall",'BCWall')
t = Internal.addGhostCells(t,t,2,adaptBCs=1,fillCorner=1)
test.testT(t,2)
# geometrical extrapolation of corner cells
t = C.newPyTree(['Base',a,b,c,d])
t = C.initVars(t, '{F}=3*{CoordinateX}+2*{CoordinateY}')
t = C.initVars(t, '{centers:G}={centers:CoordinateZ}')
t = X.connectMatch(t,dim=3)
t = C.fillEmptyBCWith(t,"wall",'BCWall')
t = Internal.addGhostCells(t,t,2,adaptBCs=0,fillCorner=0)
test.testT(t,3)
# geometrical extrapolation of corner cells
t = C.newPyTree(['Base',a,b,c,d])
t = C.initVars(t, '{F}=3*{CoordinateX}+2*{CoordinateY}')
t = C.initVars(t, '{centers:G}={centers:CoordinateZ}')
t = X.connectMatch(t,dim=3)
t = C.fillEmptyBCWith(t,"wall",'BCWall')
t = Internal.addGhostCells(t,t,2,adaptBCs=1,fillCorner=0)
test.testT(t,4)
