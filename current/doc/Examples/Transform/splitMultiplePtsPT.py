# - splitMultiplePts (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import Connector.PyTree as X

nk = 2
z0 = G.cart((0.,0.,0.),(0.1,0.1,1.),(10,10,nk))
z1 = T.subzone(z0,(1,1,1),(5,10,nk)); z1[0] = 'cart1'
z2 = T.subzone(z0,(5,1,1),(10,5,nk)); z2[0] = 'cart2'
z3 = T.subzone(z0,(5,5,1),(10,10,nk)); z3[0] = 'cart3'
z0 = T.translate(z0,(-0.9,0.,0.)); z0[0] = 'cart0'
z4 = G.cart((-0.9,0.9,0.),(0.1,0.1,1.),(19,5,nk)); z4[0] = 'cart4'
t = C.newPyTree(['Base',z1,z2,z3,z4])
t = X.connectMatch(t,dim=2)
t = C.fillEmptyBCWith(t, 'wall', 'BCWall', dim=2)
t = T.splitMultiplePts(t, dim=2)
C.convertPyTree2File(t, 'out.cgns')
