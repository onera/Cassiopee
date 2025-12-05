# Test des fonctions :
# - convertArrays3D2Arrays

import Generator as G
import Converter.Array3D
import KCore.test as test

a = G.cart( (0,0,0), (0.1, 0.2, 1.), (11, 4, 2))
b = Converter.Array3D.convertArrays2Arrays3D([a])
c = Converter.Array3D.convertArrays3D2Arrays(b)
test.testA(c, 1)
