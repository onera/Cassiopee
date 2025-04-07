# - computeDiff (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

def F(x,y):
    return (x-5)*(x-5)
def celln(y):
    if ( y > 5. ): return True
    else: return False

#-------------------------------
# 2D structure + raccords match
#-------------------------------
ni = 10; nj = 10; nk = 1
a = G.cart((0,0,0), (1,1,1), (ni,nj,nk))
a = C.initVars(a, 'Density', F, ['CoordinateX', 'CoordinateY'])
b = G.cart((9,0,0), (1,1,1), (ni,nj,nk))
b = C.initVars(b, 'Density', F, ['CoordinateX', 'CoordinateY'])

a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', b, 'imin',
                 trirac=[1,2])
b = C.addBC2Zone(b, 'match2', 'BCMatch', 'imin', a, 'imax',
                 trirac=[1,2])

t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t = P.computeDiff(t, 'Density')
test.testT(t,1)
