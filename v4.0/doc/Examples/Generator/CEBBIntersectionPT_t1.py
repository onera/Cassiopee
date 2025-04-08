# - CEBBIntersection (pyTree)-
import Generator.PyTree as G
import Converter.PyTree as C
import Transform.PyTree as T
import KCore.test as test

# entre deux zones
ni = 11; nj = 3; nk = 11
a1 = G.cart((0.,0.,0.), (0.1,0.1,0.2),(ni, nj,nk))
a1 = C.addBC2Zone(a1, 'wall1','BCWall','jmin')
a1 = C.addBC2Zone(a1, 'match1','BCMatch','imin',a1,'imax',[1,2,3])
a1 = C.addBC2Zone(a1, 'match2','BCMatch','imax',a1,'imin',[1,2,3])
a1 = C.fillEmptyBCWith(a1,'overlap','BCOverlap')

a2 = G.cart((0.9,0.,0.), (0.1,0.1,0.2),(ni, nj,nk))
a2 = T.rotate(a2, (0,0,0), (0,0,1), 12.)
a2 = C.addBC2Zone(a2, 'wall1','BCWall','jmin')
a2 = C.addBC2Zone(a2, 'match1','BCMatch','imin',a2,'imax',[1,2,3])
a2 = C.addBC2Zone(a2, 'match2','BCMatch','imax',a2,'imin',[1,2,3])
a2 = C.fillEmptyBCWith(a2,'overlap','BCOverlap')
res = G.CEBBIntersection(a1, a2)
test.testO(res, 1)
