# - selectCells2 (array) -
import Converter as C
import Generator as G
import Post as P
import KCore.test as test

def G(x, y, z):
    return (x + y + z > 5.)

def F(a):
    b = C.initVars(a, 'tag', G, ['x','y','z'], isVectorized=True)
    tag = C.extractVars(b, ['tag'])
    return P.selectCells2(a, tag, strict=1)

test.stdTestA(F)
