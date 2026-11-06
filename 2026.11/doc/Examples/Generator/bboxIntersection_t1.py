# - bboxIntersection (array) -
import Generator as G
import KCore.test as test

ni = 11; nj = 3; nk = 11
a1 = G.cart((0.,0.,0.), (0.1,0.1,0.2),(ni, nj,nk))
a2 = G.cart((0.5,0.05,0.01), (0.1,0.1,0.2),(ni, nj,nk))
intersect = G.bboxIntersection(a1,a2)
test.testO(intersect)
