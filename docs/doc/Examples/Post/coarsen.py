# - coarsen (array) -
import Post as P
import Converter as C
import Generator as G
import Transform as T

# coarsen all cells of a square
ni = 21; nj = 21; nk = 11
hi = 2./(ni-1); hj = 2./(nj-1); hk = 1./(nk-1)
m = G.cart((0.,0.,0.),(hi,hj,hk), (ni,nj,nk))

hi = hi/2; hj = hj/2; hk = hk/2
m2 = G.cart((0.,0.,0.),(hi,hj,hk), (ni,nj,nk))
m2 = T.subzone(m2,(3,3,6),(m[2]-2,m[3]-2,6))
m2 = T.translate(m2, (0.75,0.75,0.25))
m2 = T.perturbate(m2, 0.51)
tri = G.delaunay(m2)

npts = tri[2].shape[1]
indic = C.array('indic', npts, 1, 1)
indic = C.initVars(indic, 'indic', 1)

sol = P.coarsen(tri, indic, argqual=0.25, tol=1.e6)
C.convertArrays2File([tri, sol], 'out.plt')
