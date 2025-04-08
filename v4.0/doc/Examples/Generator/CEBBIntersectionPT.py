# - CEBBIntersection (pyTree)-
import Generator.PyTree as G
import Transform.PyTree as T

ni = 11; nj = 3; nk = 11
a1 = G.cart((0.,0.,0.), (0.1,0.1,0.2),(ni, nj,nk))
a2 = G.cart((1.,0.,0.), (0.1,0.1,0.2),(ni, nj,nk))
a2 = T.rotate(a2, (0,0,0), (0,0,1), 12.)
print(G.CEBBIntersection(a1, a2))
