# - distance2Walls (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Dist2Walls.PyTree as DTW
import Geom.PyTree as D
import KCore.test as test

a1 = G.cart((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 11)); a1[0] = 'cart1'
a2 = G.cart((1., 0., 0.), (0.1, 0.1, 0.1), (11, 21, 11)); a2[0] = 'cart2'
a3 = G.cart((0.,0.,1.), (0.1, 0.1, 0.1), (11, 21, 11)); a3[0] = 'cart3'
a2 = C.convertArray2NGon(a2)
t = C.newPyTree(['Base',a1,a2,a3])
s = D.sphere6((0,0,0), 1., N=20, ntype='TRI')
DTW._distance2Walls(t, s, type='ortho')
test.testT(t)
