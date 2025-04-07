# - computeGrad2 (pyTree) -
import Converter.Internal as CI
import Converter.PyTree   as C
import Post.PyTree        as P
import Generator.PyTree   as G
import Transform.PyTree   as T
import Connector.PyTree   as X
import KCore.test         as test

# Cas 3D NGon avec raccord match
# ===============================
ni = 30
nj = 30
nk = 10
lx = 1.
ly = 1.
lz = 0.25
dx = lx/float(ni-1)
dy = ly/float(nj-1)
dz = lz/float(nk-1)

a = G.cartNGon((1,1,1),           (dx,dy,dz), (ni,nj,nk)); a[0]='cart1'
b = G.cartNGon((lx+1,1+(ni/2.)*dy,1), (dx,dy,dz), (ni,nj,nk)); b[0]='cart2'
c = G.cartNGon((lx+1,1-(ni/2.-1)*dy,1), (dx,dy,dz), (ni,nj,nk)); c[0]='cart3'

t = C.newPyTree(['Base',a,b,c])

t = C.initVars(t, '{centers:F}={centers:CoordinateX}*{centers:CoordinateY}*{centers:CoordinateZ}')
t = C.initVars(t, '{G}=2.3')
t = C.initVars(t, '{centers:H}={centers:CoordinateY}')
t = C.initVars(t, '{centers:M}={centers:CoordinateX}')
t = X.connectMatch(t)
t = C.fillEmptyBCWith(t,"wall",'BCWall')

P._computeGrad2(t, 'centers:F')

test.testT(t,1)
