# - selectCells2 (array) -
import Converter as C
import Generator as G
import Post as P
import KCore.test as test

def G(x, y, z):
    if (x+y+z > 5.): return True
    else: return False

def F(a):
    b = C.initVars(a, 'tag', G, ['x','y','z'])
    tag = C.extractVars(b, ['tag'])
    return P.selectCells2(a, tag, strict=1)

test.stdTestA(F)
