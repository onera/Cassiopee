# - checkPointInCEBB (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import KCore.test as test

#------------------
# CAS 3D STRUCTURE
#------------------
Ni = 20; Nj = 20
a2 = G.cart((-0.1,0,0),(0.5/Ni, 0.5/Nj, 1), (Ni,Nj,2))
a2 = T.rotate(a2, (-0.1,0,0), (0,0,1), 0.22)
a2 = C.addBC2Zone(a2, 'wall1','BCWall','jmin')
a2 = C.addBC2Zone(a2, 'match1','BCMatch','imin',a2,'imax',[1,2,3])
a2 = C.addBC2Zone(a2, 'match2','BCMatch','imax',a2,'imin',[1,2,3])
a2 = C.fillEmptyBCWith(a2,'overlap','BCOverlap')
val = G.checkPointInCEBB(a2, (0.04839, 0.03873, 0.5))
test.testO(val,1)

#------------------
# CAS 2D STRUCTURE
#------------------
Ni = 20; Nj = 20
a2 = G.cart((-0.1,0,0),(0.5/Ni, 0.5/Nj, 1), (Ni,Nj,1))
a2 = T.rotate(a2, (-0.1,0,0), (0,0,1), 0.22)
a2 = C.addBC2Zone(a2, 'wall1','BCWall','jmin')
a2 = C.addBC2Zone(a2, 'match1','BCMatch','imin',a2,'imax',[1,2,3])
a2 = C.addBC2Zone(a2, 'match2','BCMatch','imax',a2,'imin',[1,2,3])
a2 = C.fillEmptyBCWith(a2,'overlap','BCOverlap')
val = G.checkPointInCEBB(a2, (0.04839, 0.03873, 0.))
test.testO(val,2)
