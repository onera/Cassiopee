# - getFamilyBCNamesOfType (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0.,0.,0), (0.01,0.01,1.), (20,20,2))
b = G.cart((1.,0.,0), (0.01,0.01,1.), (20,20,2))

a = C.addBC2Zone(a, 'walla', 'FamilySpecified:CARTER', 'imin')
a = C.addBC2Zone(a, 'nref', 'FamilySpecified:LOIN', 'imax')
b = C.addBC2Zone(b, 'wallb', 'FamilySpecified:CARTER', 'jmin')

t = C.newPyTree(['Base',a,b])

t[2][1] = C.addFamily2Base(t[2][1], 'CARTER', bndType='BCWall')
t[2][1] = C.addFamily2Base(t[2][1], 'LOIN', bndType='BCFarfield')

# Toutes les familyBCs de type BCwall
names1 = C.getFamilyBCNamesOfType(t, 'BCWall')
# Toutes les familyBCs de l'arbre
names2 = C.getFamilyBCNamesOfType(t)
names1.sort(); names2.sort()
test.testO([names1, names2], 1)
