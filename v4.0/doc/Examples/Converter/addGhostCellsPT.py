# - addGhostCells (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

ni = 11; nj = 11; nk = 8
a = G.cart((0,0,0), (1.,1.,1.), (ni,nj,nk)); a[0]='cart1'
b = G.cart((0,0,-3.5), (1.,1.,0.5), (ni,nj,nk)); b[0]='cart2'
a = C.addBC2Zone(a,'match','BCMatch','kmin',b[0],[1,ni,1,nj,nk,nk],[1,2,3])
b = C.addBC2Zone(b,'match','BCMatch','kmax',a[0],[1,ni,1,nj,1,1],[1,2,3])
t = C.newPyTree(['Base',a,b])
t = C.addBC2Zone(t, 'wall', 'BCWall', 'imin')
t = C.initVars(t, '{F}=3*{CoordinateX}+2*{CoordinateY}')

t = C.addGhostCells(t, t, 2, adaptBCs=1)
C.convertPyTree2File(t, 'out.cgns')
