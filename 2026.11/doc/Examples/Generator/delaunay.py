# - delaunay (array) -
import Generator as G
import Converter as C
ni = 11; nj = 11; nk = 1
hi = 1./(ni-1); hj = 1./(nj-1); hk = 1.
a = G.cart((0.,0.,0.), (hi,hj,hk), (ni,nj,nk))
b = G.delaunay(a)
C.convertArrays2File([a,b], "out.plt")
