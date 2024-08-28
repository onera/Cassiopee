# - dist2wallIBM (pyTree) -
import Generator.PyTree as G
import Geom.IBM as DIBM
import Connector.IBM as XIBM
import Converter.PyTree as C
import Converter.Internal as Internal
import Geom.PyTree as D

N = 50; h = 1./(N-1)
tb = D.sphere((0.,0.,0.),0.2,2*N); tb[0] = 'sphere'
DIBM._setIBCType(tb, 'Musker')

t = G.cart((-1,-1,-1), (h,h,h), (2*N,2*N,2*N)); t[0] = 'cart'
C._fillEmptyBCWith(t, 'farfield', 'BCFarfield')
XIBM._dist2wallIBM(t, tb, dimPb=3)
C.convertPyTree2File(t, 'out.cgns')
