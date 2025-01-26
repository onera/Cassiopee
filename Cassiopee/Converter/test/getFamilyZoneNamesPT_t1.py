# - getFamilyZoneNames (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0.,0.,0), (0.01,0.01,1.), (20,20,2))
b = G.cart((1.,0.,0), (0.01,0.01,1.), (20,20,2))
C._tagWithFamily(a, 'CARTER')
C._tagWithFamily(b, 'CARTER')

t = C.newPyTree(['Base',a,b])
C._addFamily2Base(t[2][1], 'CARTER')

# Toutes les family zone names de l'arbre
names = C.getFamilyZoneNames(t)
test.testO(names, 1)
