# - coarsen (pyTree)-
import Post.PyTree as P
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T

ni = 21; nj = 21; nk = 1
hi = 2./(ni-1); hj = 2./(nj-1)
m = G.cart((0.,0.,0.),(hi,hj,1.), (ni,nj,nk)); m = T.perturbate(m, 0.51)
tri = G.delaunay(m); tri = C.initVars(tri, 'centers:indic', 1.)

sol = P.coarsen(tri, 'indic'); sol[0] = 'coarse'

C.convertPyTree2File([sol,tri], 'out.cgns')
