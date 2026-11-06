# - dist2wallIBM (pyTree) -
import Generator.PyTree as G
import Geom.IBM as DIBM
import Connector.IBM as XIBM
import Converter.PyTree as C
import Converter.Internal as Internal
import KCore.test as test

# --- Front Type=1
N = 50; h = 1./(N-1)
tb = G.cylinder((0.,0.,0.), 0.2, 1., 360., 0., 10., (2*N,1,1)); tb[0] = 'cyl'
DIBM._setIBCType(tb, 'Musker')

t = G.cart((-1,-1,0), (h,h,1), (2*N,2*N,1)); t[0] = 'cart'
C._fillEmptyBCWith(t, 'farfield', 'BCFarfield')
XIBM._dist2wallIBM(t, tb, dimPb=2)
test.testT(t, 1)

# --- Front Type = 42
N = 50; h = 1./(N-1)
tb = G.cylinder((0.,0.,0.), 0.2, 1., 360., 0., 10., (2*N,1,1)); tb[0] = 'cyl'
DIBM._setIBCType(tb, 'Musker')

t = G.cart((-1,-1,0), (h,h,1), (2*N,2*N,1)); t[0] = 'cart'
C._fillEmptyBCWith(t, 'farfield', 'BCFarfield')
XIBM._dist2wallIBM(t, tb, dimPb=2)
test.testT(t, 2)
