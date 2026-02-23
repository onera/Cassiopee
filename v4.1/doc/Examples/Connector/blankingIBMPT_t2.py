# - blankingIBM 3D (pyTree) -
import Generator.PyTree as G
import Geom.IBM as DIBM
import Connector.IBM as XIBM
import Converter.PyTree as C
import Converter.Internal as Internal
import Geom.PyTree as D
import KCore.test as test

N = 30; h = 1./(N-1)
s1 = D.sphere((0.,0.,0.),0.2,N//2); s1[0] = 'sphere'
DIBM._setIBCType(s1, 'Musker')
tb = C.newPyTree(['SPH', Internal.getZones(s1)])

a = G.cart((-1,-1,-1), (2*h,2*h,2*h), (N,N,N)); a[0] = 'cart'
t = C.newPyTree(['CART', Internal.getZones(a)])
C._fillEmptyBCWith(t, 'farfield', 'BCFarfield')

XIBM._dist2wallIBM(t, tb, dimPb=3)
XIBM._blankingIBM(t, tb, dimPb=3)
test.testT(t, 1)

# --- Front Type=42
N = 30; h = 1./(N-1)
s1 = D.sphere((-0.5,-0.5,-0.5),0.2,N//2); s1[0] = 'sphere'
tb = C.newPyTree(['SPH', Internal.getZones(s1)])
DIBM._setIBCType(tb, 'noslip')

a = G.cart((-1,-1,-1), (2*h,2*h,2*h), (N,N,N)); a[0] = 'cart'
t = C.newPyTree(['CART', Internal.getZones(a)])
C._fillEmptyBCWith(t, 'farfield', 'BCFarfield')

XIBM._dist2wallIBM(t, tb, dimPb=3, frontType=42)
XIBM._blankingIBM(t, tb, dimPb=3, frontType=42)
test.testT(t, 2)

# --- Front Type = 42 , correction multicorps = True
N = 30; h = 1./(N-1)
s1 = D.sphere((-0.5,-0.5,-0.5),0.2,N//2)
s2 = D.sphere((0.5,0.5,0.5),0.2,N//2)
tb = C.newPyTree(['SPH1', Internal.getZones(s1),
                  'SPH2', Internal.getZones(s2)])
DIBM._setIBCType(tb, 'Musker')

a = G.cart((-1,-1,-1), (2*h,2*h,2*h), (N,N,N)); a[0] = 'cart'
t = C.newPyTree(['CART', Internal.getZones(a)])
C._fillEmptyBCWith(t, 'farfield', 'BCFarfield')

XIBM._dist2wallIBM(t, tb, dimPb=3, correctionMultiCorpsF42=True, frontType=42)
XIBM._blankingIBM(t, tb, dimPb=3, correctionMultiCorpsF42=True, frontType=42)
test.testT(t, 3)

# --- Front Type=42 , blankingF42=True
N = 30; h = 1./(N-1)
s1 = D.sphere((-0.5,-0.5,-0.5),0.2,N//2); s1[0] = 'sphere'
tb = C.newPyTree(['SPH', Internal.getZones(s1)])
DIBM._setIBCType(tb, 'Musker')

a = G.cart((-1,-1,-1), (2*h,2*h,2*h), (N,N,N)); a[0] = 'cart'
t = C.newPyTree(['CART', Internal.getZones(a)])
C._fillEmptyBCWith(t, 'farfield', 'BCFarfield')

XIBM._dist2wallIBM(t, tb, dimPb=3, frontType=42)
XIBM._blankingIBM(t, tb, dimPb=3, frontType=42, blankingF42=True)
test.testT(t, 4)

# --- Front Type=42 , twoFronts=True
N = 30; h = 1./(N-1)
s1 = D.sphere((-0.5,-0.5,-0.5),0.2,N//2); s1[0] = 'sphere'
tb = C.newPyTree(['SPH', Internal.getZones(s1)])
DIBM._setIBCType(tb, 'Musker')

a = G.cart((-1,-1,-1), (2*h,2*h,2*h), (N,N,N)); a[0] = 'cart'
t = C.newPyTree(['CART', Internal.getZones(a)])
C._fillEmptyBCWith(t, 'farfield', 'BCFarfield')

XIBM._dist2wallIBM(t, tb, dimPb=3, frontType=42)
XIBM._blankingIBM(t, tb, dimPb=3, frontType=42, twoFronts=True)
test.testT(t, 5)
