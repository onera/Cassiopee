# - coarsen (pyTree)-
import Post.PyTree as P
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

# test sur une zone
ni = 21; nj = 21; nk = 1
hi = 2./(ni-1); hj = 2./(nj-1)
m = G.cart((0.,0.,0.),(hi,hj,1.), (ni,nj,nk)); m = T.perturbate(m, 0.51)
tri = G.delaunay(m); tri = C.initVars(tri, 'centers:indic', 1.)
tri = C.initVars(tri, 'G', 2.)
sol = P.coarsen(tri, 'indic')
test.testT(sol, 1)

# test sur une base
t = C.newPyTree(['Base',2]); t[2][1][2] += [tri]
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
sol = P.coarsen(t[2][1], 'indic')
test.testT(sol, 2)

# test sur un arbre
sol = P.coarsen(t, 'indic')
test.testT(sol,3)
