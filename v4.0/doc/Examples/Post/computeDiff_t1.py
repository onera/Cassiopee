# - computeDiff (array) -
import Converter as C
import Post as P
import Generator as G
import KCore.test as test

def F(x):
    if ( x > 5. ): return True
    else : return False

def celln(y):
    if ( y > 5. ): return True
    else : return False

#--------------
# 2D structure
#--------------
ni = 30; nj = 40; nk = 1
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m, 'ro', F, ['x'])
p = P.computeDiff(m, 'ro')
test.testA([p],1)

m = C.initVars(m, 'cellN', celln, ['y'])
p = P.computeDiff(m,'ro')
test.testA([p],2)

#--------------
# 3D structure
#--------------
ni = 30; nj = 40; nk = 11
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),10./(nk-1)), (ni,nj,nk))
m = C.initVars(m, 'ro', F, ['x'])
p = P.computeDiff(m,'ro')
test.testA([p],3)
#
m = C.initVars(m, 'cellN', celln, ['y'])
p = P.computeDiff(m,'ro')
test.testA([p],4)
