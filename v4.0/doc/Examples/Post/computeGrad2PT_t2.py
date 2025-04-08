# - computeGrad2 (pyTree) -
import Converter.PyTree   as C
import Post.PyTree        as P
import Generator.PyTree   as G
import Transform.PyTree   as T
import Connector.PyTree   as X
import KCore.test         as test

# Cas 2D avec raccord match
# =========================
ni = 20
nj = 20
nk =  1
dx = 1./float(ni-1)
dy = 1./float(nj-1)
dz = 1.

a = G.cart((1,1,1), (dx,dy,dz), (ni,nj,nk)); a[0]='cart1'
b = G.cart((2,1.-5*dx,1), (dx,dy,dz), (ni,nj+10,nk)); b[0]='cart2'

a = T.reorder(a,(-2,1,3))

a = C.initVars(a, '{centers:var}=103.')
b = C.initVars(b, '{centers:var}=104.')

t = C.newPyTree(['Base',a,b])
t = X.connectMatch(t,dim=2)
t = C.initVars(t, '{centers:F}={centers:CoordinateX}*{centers:CoordinateY}')

P._computeGrad2(t, 'centers:F')

test.testT(t,1)
