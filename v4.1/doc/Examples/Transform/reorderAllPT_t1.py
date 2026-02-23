# - reorderAll (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

def F(x,y): return x*x + 2*y
ni = 30; nj = 40
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
a = C.initVars(a, 'centers:cellN', 1.)
a = C.initVars(a, 'Density', F, ['CoordinateX','CoordinateY'])
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imax')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmax')

#-------------------------------------------------
# test 1 : blocs recouvrants pareillement orientes
#-------------------------------------------------
b = T.rotate(a, (0.2,0.2,0.), (0.,0.,1.), 15.); b[0] = 'cart2'
B = T.reorderAll([a,b], 1)
test.testT(B, 1)

#-------------------------------------------------
# test 2 : blocs recouvrants orientes differemment
#-------------------------------------------------
b = T.reorder(b, (-1,2,3))
B = T.reorderAll([a,b], 1)
test.testT(B, 2)

#---------------------------------
# test 3: blocs sans intersection
#---------------------------------
b = T.translate(a, (0.,12,0.))
A = [a,b]
B = T.reorderAll(A)
test.testT(B, 3)
