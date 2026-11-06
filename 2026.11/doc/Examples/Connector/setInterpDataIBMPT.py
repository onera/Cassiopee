# - setInterpDataIBM (pyTree) -
import Generator.PyTree as G
import Geom.IBM as DIBM
import Connector.IBM as XIBM
import Converter.PyTree as C
import Converter.Internal as Internal
import Geom.PyTree as D
import Transform.PyTree as T

N = 50; h = 1./(N-1)
s = D.sphere((0.,0.,0.),0.2,2*N); s[0] = 'sphere'
DIBM._setIBCType(s, 'Musker')
tb = C.newPyTree(['SPH', Internal.getZones(s)])
C._addState(tb, 'GoverningEquations', 'NSTurbulent')

a = G.cart((-1,-1,-1), (h,h,h), (2*N,2*N,2*N)); a[0] = 'cart'
t = C.newPyTree(['CART', Internal.getZones(a)])
C._fillEmptyBCWith(t, 'farfield', 'BCFarfield')

XIBM._dist2wallIBM(t, tb, dimPb=3)
XIBM._blankingIBM(t, tb, dimPb=3)
tc = C.node2Center(t)
t, tc, front, front2 = XIBM.buildFrontIBM(t, tc, dimPb=3, frontType=42, check=True)
XIBM._setInterpDataIBM(t, tc, tb, front, frontType=42, dimPb=3)

C.convertPyTree2File(t, 'out.cgns')
