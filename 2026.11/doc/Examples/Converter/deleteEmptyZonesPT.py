# - deleteEmptyZones (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P

a = G.cart((0,0,0), (1,1,1), (3,3,3))
b = P.selectCells(a, '{CoordinateX} > 12')
c = P.selectCells(a, '{CoordinateX} > 15')

t = C.newPyTree(['Base',c,a,b])
C._deleteEmptyZones(t)
C.convertPyTree2File(t, 'out.cgns')
