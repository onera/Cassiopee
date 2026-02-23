# - chimeraInfo (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X

a = G.cylinder((0,0,0),1.,3.,360,0,1,(200,30,4)); a[0] = 'cylindre1'
C._addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
C._addBC2Zone(a, 'ov1', 'BCOverlap', 'jmax')
b = G.cylinder((4,0,0),1.,3.,360,0,1,(200,30,4)); b[0] = 'cylindre2'
C._addBC2Zone(b, 'wall1', 'BCWall', 'jmin')
C._addBC2Zone(b, 'ov1', 'BCOverlap', 'jmax')
t = C.newPyTree(['Cyl1',a,'Cyl2',b])
t = X.connectMatch(t, dim=3)
C._fillEmptyBCWith(t,'nref','BCFarfield', dim=3)
C._initVars(t[2][2],'centers:cellN',2)
t = X.setInterpolations(t, loc='cell', double_wall=1, storage='direct')
X._chimeraInfo(t, type='orphan')
orphanPts = X.extractChimeraInfo(t,type='orphan')
C.convertPyTree2File(orphanPts, "orphanPts.cgns")
X._chimeraInfo(t,type='extrapolated')
out = X.extractChimeraInfo(t,type='cf>1.5')
C.convertPyTree2File(out,"out.cgns")
