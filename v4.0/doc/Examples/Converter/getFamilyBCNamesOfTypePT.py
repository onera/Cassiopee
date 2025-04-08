# - getFamilyBCNamesOfType (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0.,0.,0), (0.01,0.01,1.), (20,20,2))
b = G.cart((1.,0.,0), (0.01,0.01,1.), (20,20,2))

a = C.addBC2Zone(a, 'walla', 'FamilySpecified:CARTER', 'imin')
b = C.addBC2Zone(b, 'wallb', 'FamilySpecified:CARTER', 'jmin')

t = C.newPyTree(['Base',a,b])

C._addFamily2Base(t[2][1], 'CARTER', bndType='BCWall')

# Toutes les familyBCs de type BCwall
names = C.getFamilyBCNamesOfType(t, 'BCWall'); print(names)
#>> ['CARTER']
# Toutes les familyBCs de l'arbre
names = C.getFamilyBCNamesOfType(t); print(names)
#>> ['CARTER']
