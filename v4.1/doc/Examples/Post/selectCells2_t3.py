# - selectCells2 (array) -
import Converter as C
import Generator as G
import Post as P
import KCore.test as test

def F(x, y, z):
    if (x+2*y+z > 20.): return True
    else: return False

# List
a    = G.cart( (0,0,0), (1,1,1), (11,11,11) )
taga = C.initVars(a, 'tag', F, ['x','y','z'])
taga = C.extractVars(taga,['tag'])

b    = G.cart( (1,1,1), (1,1,1), (11,11,11) )
tagb = C.initVars(b, 'tag', F, ['x','y','z'])
tagb = C.extractVars(tagb,['tag'])

c    = [a,b]
tag  = [taga,tagb]

c = P.selectCells2(c, tag)

test.testA([a], 1)
