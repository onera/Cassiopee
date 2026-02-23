# - extractBCOfName (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cylinder((0,0,0), 1., 1.5, 360., 0., 1., (100,30,10))
C._initVars(a, '{F1}={CoordinateX}');C._initVars(a, '{centers:G1}={centers:CoordinateY}')
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'walla', 'FamilySpecified:CARTER', 'imin')
t = C.newPyTree(['Base',3,'Skin',2]); t[2][1][2] += [a]
t[2][1] = C.addFamily2Base(t[2][1], 'CARTER', bndType='BCWall')

Z1 = C.extractBCOfName(a, 'wall1')
Z2 = C.extractBCOfName(a, 'FamilySpecified:CARTER')
t[2][2][2] += Z1+Z2
test.testT(t,1)

# Version wild card
a = G.cylinder((0,0,0), 1., 1.5, 360., 0., 1., (100,30,10))
C._initVars(a, '{F1}={CoordinateX}');C._initVars(a, '{centers:G1}={centers:CoordinateY}')
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'walla', 'FamilySpecified:CARTER', 'imin')
t = C.newPyTree(['Base',3,'Skin',2]); t[2][1][2] += [a]
t[2][1] = C.addFamily2Base(t[2][1], 'CARTER', bndType='BCWall')
zones = C.extractBCOfName(t, 'wall*')
test.testT(zones,2)
