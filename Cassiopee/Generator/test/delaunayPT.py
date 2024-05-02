# - delaunay (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

ni = 11; nj = 11; nk = 1
hi = 1./(ni-1); hj = 1./(nj-1); hk = 1.
a = G.cart((0.,0.,0.), (hi,hj,hk), (ni,nj,nk))
b = G.delaunay(a); b[0] = 'delaunay'
t = C.newPyTree(['Base', 2, a,b])
C.convertPyTree2File(t, 'out.cgns')
