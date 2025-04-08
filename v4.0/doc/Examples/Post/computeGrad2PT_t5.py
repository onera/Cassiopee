# - computeGrad2 (pyTree) -
import Converter.PyTree   as C
import Post.PyTree        as P
import Generator.PyTree   as G
import KCore.test         as test

# Cas 3D NGon avec volume calcule
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

a = G.cartNGon((1,1,1),(dx,dy,dz),(ni,nj,nk)); a[0]='cart1'

t = C.newPyTree(['Base',a])
G._getVolumeMap(t)

t = C.initVars(t, '{centers:F}={centers:CoordinateX}*{centers:CoordinateY}*{centers:CoordinateZ}')

P._computeGrad2(t, 'centers:F')

test.testT(t,1)
