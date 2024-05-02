# - extractBCOfName (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cylinder((0,0,0), 1., 1.5, 360., 0., 1., (100,30,10))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'walla', 'FamilySpecified:CARTER', 'imin')
t = C.newPyTree(['Base',3,'Skin',2]); t[2][1][2] += [a]
t[2][1] = C.addFamily2Base(t[2][1], 'CARTER', bndType='BCWall')

Z1 = C.extractBCOfName(a, 'wall*')
Z2 = C.extractBCOfName(a, 'FamilySpecified:CARTER')
C.convertPyTree2File(Z1+Z2, 'out.cgns')
