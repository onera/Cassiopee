# - CEBBIntersection (array) -
import Generator as G
import Converter as C
import Transform as T
import KCore.test as test

ni = 11; nj = 3; nk = 11
a1 = G.cart((0.,0.,0.), (0.1,0.1,0.2),(ni, nj,nk))
a2 = G.cart((1.,0.,0.), (0.1,0.1,0.2),(ni, nj,nk))
a2 = T.rotate(a2, (0,0,0), (0,0,1), 12.)
res = G.CEBBIntersection(a1, a2)
test.testO(res, 1)

a2 = T.translate(a2, (-0.1,-0.1,0))
res = G.CEBBIntersection(a1, a2)
test.testO(res, 2)
