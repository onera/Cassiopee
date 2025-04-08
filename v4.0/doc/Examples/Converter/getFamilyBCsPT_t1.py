# - getFamilyBCs (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0.,0.,0), (0.01,0.01,1.), (20,20,2))
b = G.cart((1.,0.,0), (0.01,0.01,1.), (20,20,2))

a = C.addBC2Zone(a, 'walla', 'FamilySpecified:CARTER', 'imin')
b = C.addBC2Zone(b, 'wallb', 'FamilySpecified:CARTER', 'jmin')

t = C.newPyTree(['Base']); t[2][1][2] += [a,b]

t[2][1] = C.addFamily2Base(t[2][1], 'CARTER', bndType='BCWall')

B1 = C.getFamilyBCs(t, 'CARTER')
test.testO(B1, 1)
