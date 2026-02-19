# - fillEmptyBC (pyTree) -

# NGON
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# NGON3
a = G.cartNGon((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))
# ajoute en tant que face list
a = C.fillEmptyBCWith(a, 'wallv', 'BCWallViscous')
t = C.newPyTree(['Base', a])
test.testT(t, 1)

# NGON4
a = G.cartNGon((0.,0.,0.),(0.1,0.1,0.1),(10,10,10), api=3)
# ajoute en tant que face list
a = C.fillEmptyBCWith(a, 'wallv', 'BCWallViscous')
t = C.newPyTree(['Base', a])
test.testT(t, 2)
