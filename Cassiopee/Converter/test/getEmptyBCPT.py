# - getEmptyBC (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

a1 = G.cart((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 2)); a1[0] = 'cart1'
a1 = C.addBC2Zone(a1, 'wall1', 'BCWall', 'imin')
a2 = G.cart((1., 0.2, 0.), (0.1, 0.1, 0.1), (11, 21, 2)); a2[0] = 'cart2'
a2 = C.addBC2Zone(a2, 'wall1', 'BCWall', 'imax')
t = C.newPyTree(['Base',a1,a2])
# Returns undefined windows (as range list since structured)
wins = C.getEmptyBC(t,2); print(wins)
#>> [[ [[11, 11, 1, 21, 1, 2], ..., [1, 11, 21, 21, 1, 2]], [[1, 1, 1, 21, 1, 2], ..., [1, 11, 21, 21, 1, 2]] ]]
# Returns undefined windows (as face indices list)
t = C.convertArray2NGon(t)
faceList = C.getEmptyBC(t,2); print(faceList)
#>> [[ [array([ 11, 230, ...], dtype=int32)], [array([  1, 221, 222, 632, ...], dtype=int32)] ]]
C.convertPyTree2File(t, 'out.cgns')
