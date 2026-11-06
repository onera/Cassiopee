# - bboxIntersection (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# 3D STRUCTURE
ni = 11; nj = 3; nk = 11
a1 = G.cart((0.,0.,0.), (0.1,0.1,0.2),(ni, nj,nk))
a2 = G.cart((0.5,0.05,0.01), (0.1,0.1,0.2),(ni, nj,nk))
a1 = C.addBC2Zone(a1, 'wall1','BCWall','jmin')
a1 = C.addBC2Zone(a1, 'match1','BCMatch','imin',a1,'imax',[1,2,3])
a1 = C.addBC2Zone(a1, 'match2','BCMatch','imax',a1,'imin',[1,2,3])
a1 = C.fillEmptyBCWith(a1,'overlap','BCOverlap')
a2 = C.addBC2Zone(a2, 'wall1','BCWall','jmin')
a2 = C.addBC2Zone(a2, 'match1','BCMatch','imin',a2,'imax',[1,2,3])
a2 = C.addBC2Zone(a2, 'match2','BCMatch','imax',a2,'imin',[1,2,3])
a2 = C.fillEmptyBCWith(a2,'overlap','BCOverlap')
intersect = G.bboxIntersection(a1, a2)
test.testO(intersect,1)

# TETRA
at1 = C.convertArray2Tetra(a1)
intersect = G.bboxIntersection(a1, a2)
test.testO(intersect,2)

# HEXA
at1 = C.convertArray2Hexa(a1)
intersect = G.bboxIntersection(a1, a2)
test.testO(intersect,3)

# 2D STRUCTURE
ni = 11; nj = 3; nk = 1
a1 = G.cart((0.,0.,0.), (0.1,0.1,1.),(ni, nj,nk))
a2 = G.cart((0.5,0.05,0.), (0.1,0.1,1.),(ni, nj,nk))
a1 = C.addBC2Zone(a1, 'wall1','BCWall','jmin')
a1 = C.addBC2Zone(a1, 'match1','BCMatch','imin',a1,'imax',[1,2,3])
a1 = C.addBC2Zone(a1, 'match2','BCMatch','imax',a1,'imin',[1,2,3])
a2 = C.addBC2Zone(a2, 'wall1','BCWall','jmin')
a2 = C.addBC2Zone(a2, 'match1','BCMatch','imin',a2,'imax',[1,2,3])
a2 = C.addBC2Zone(a2, 'match2','BCMatch','imax',a2,'imin',[1,2,3])
intersect = G.bboxIntersection(a1, a2)
test.testO(intersect,4)

# TRI
at1 = C.convertArray2Tetra(a1)
intersect = G.bboxIntersection(a1, a2)
test.testO(intersect,5)

# QUAD
at1 = C.convertArray2Hexa(a1)
intersect = G.bboxIntersection(a1, a2)
test.testO(intersect,6)
