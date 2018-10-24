# - computeGrad (array) -
import Converter as C
import Post as P
import Generator as G
import KCore.test as T

def F(x,y): return 2*x+x*y

# cas 3D structure
ni = 30; nj = 40; nk = 3
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m, 'ro', F, ['x','y'])
p = P.computeGrad(m, 'ro') # p is defined on centers
T.testA([p], 1)

# cas 2D structure
ni = 10; nj = 20; nk = 1
m2 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m2 = C.initVars(m2, 'ro', F, ['x','y'])
p2 = P.computeGrad(m2, 'ro') # p is defined on centers
T.testA([p2], 2)

# test sur une liste
p3 = P.computeGrad([m,m2], 'ro') # p is defined on centers
T.testA(p3, 3)
