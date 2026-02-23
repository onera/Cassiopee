# - rmGhostCells (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal

a = G.cart((1,1,1), (1.,1.,1.), (4,2,3)); a[0]='cart1'
b = G.cart((1,1,-3), (1.,1.,0.5), (4,2,9)); b[0]='cart2'
a = C.addBC2Zone(a,'match','BCMatch','kmin',b[0],[1,4,1,2,9,9],[1,2,3])
b = C.addBC2Zone(b,'match','BCMatch','kmax',a[0],[1,4,1,2,1,1],[1,2,3])
t = C.newPyTree(['Base',a,b])

t = C.addBC2Zone(t, 'wall', 'BCWall', 'imin')
t = C.initVars(t, '{F}=3*{CoordinateX}+2*{CoordinateY}')

a = t[2][1][2][0]
ag = Internal.addGhostCells(t, a, 2, adaptBCs=1)
t[2][1][2][0] = ag
ag = C.rmGhostCells(t, ag, 2, adaptBCs=1)
t[2][1][2][0] = ag
C.convertPyTree2File(t,'out.cgns')
