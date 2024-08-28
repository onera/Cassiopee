# - selectCells (array) -
import Converter as C
import Generator as G
import Post as P

a = G.cart( (0,0,0), (1,1,1), (11,11,11) )
def F(x, y, z):
    if (x + 2*y + z > 20.): return True
    else: return False

b = P.selectCells(a, F, ['x', 'y', 'z'])
c = P.selectCells(a, F, ['x', 'y', 'z'], strict=1)
d = P.selectCells(a, '{x}+2*{y}+{z}>20')
e = P.selectCells(a, '({x}>2) & ({y}>2)')
C.convertArrays2File([b,c,d,e], 'out.plt')
