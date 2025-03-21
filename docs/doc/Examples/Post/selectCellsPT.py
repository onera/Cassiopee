# - selectCells (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P

def F(x, y, z):
    if (x + 2*y + z > 20.): return True
    else: return False

a = G.cart( (0,0,0), (1,1,1), (11,11,11) )
a = P.selectCells(a, F, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
C.convertPyTree2File(a, 'out.cgns')
