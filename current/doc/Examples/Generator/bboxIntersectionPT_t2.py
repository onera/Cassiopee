# - bboxIntersection (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

boxA = G.cartTetra((2,-1,0.5), (1,2,0.5), (2,2,2))

boxB = G.cartTetra((2,-1,0.5), (2,2,1), (2,2,2))
boxB  = T.rotate(boxB,(0,0,0),(10,5,-20))
AABB = G.BB(boxA)
OBB = G.BB(boxB, method='OBB')

intersect = G.bboxIntersection(AABB, AABB, tol=1e-10,isBB=True, method='AABB')
test.testO(intersect,1)

intersect = G.bboxIntersection(AABB, OBB, tol=1e-10, isBB=True, method='AABBOBB')
test.testO(intersect,2)

intersect = G.bboxIntersection(AABB, OBB, tol=1e-10, isBB=True, method='OBB')
test.testO(intersect,3)
