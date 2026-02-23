# test - enforceCurvature2 (pyTree)
import Converter.PyTree as C
import Geom.PyTree as D
import Generator.PyTree as G
import KCore.test as test

# Cas 1 : courbure constante
a = D.circle((0,0,0),0.1,N=51)
ni = 101; db = G.cart((0,0,0),(1./(ni-1),1,1),(ni,1,1))
db = C.addBC2Zone(db, 'wall1', 'BCWall', 'jmin')
db = C.addBC2Zone(db, 'match1', 'BCMatch', 'imin',db,'imax',[1,2])
db = C.addBC2Zone(db, 'match2', 'BCMatch', 'imax',db,'imin',[1,2])
db = C.addBC2Zone(db, 'wall2','BCWall','jmax')
db = C.addVars(db,'Density'); db = C.addVars(db,'centers:cellN')
db2 = G.enforceCurvature2(db, a)
test.testT([db2],1)

# Cas 2 : courbure variable
a = D.naca(12.)
ni = 101; db = G.cart((0,0,0),(1./(ni-1),1,1),(ni,1,1))
db = C.addBC2Zone(db, 'wall1','BCWall','jmin')
db = C.addBC2Zone(db, 'match1','BCMatch','imin',db,'imax',[1,2])
db = C.addBC2Zone(db, 'match2','BCMatch','imax',db,'imin',[1,2])
db = C.addBC2Zone(db, 'wall2','BCWall','jmax')
db = C.addVars(db,'Density'); db = C.addVars(db,'centers:cellN')
db2 = G.enforceCurvature2(db, a)
test.testT([db2],2)
