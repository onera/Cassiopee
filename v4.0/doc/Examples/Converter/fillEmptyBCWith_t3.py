# - fillEmptyBC (pyTree) -
# Non structure Element
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cartHexa((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))
# Pour l'instant ajoute en tant que face list
a = C.fillEmptyBCWith(a, 'wallv', 'BCWallViscous')
t = C.newPyTree(['Base',a])
test.testT(t)
