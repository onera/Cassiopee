# - selectCells (array) -
import Generator as G
import Post as P
import KCore.test as test

def F(x, y, z):
    if (x+2*y+z > 20.): return True
    else: return False

#
a = G.cart( (0,0,0), (1,1,1), (11,11,11) )
a = P.selectCells(a, F, varStrings=['x','y','z'] )
test.testA([a], 1)
