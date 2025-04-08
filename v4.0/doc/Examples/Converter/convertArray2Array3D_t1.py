# Test des fonctions :
# - convertArrays2Arrays3D

import Generator as G
import Converter.Array3D
import KCore.test as T

a = G.cart( (0,0,0), (0.1, 0.2, 1.), (11, 4, 1))
b = Converter.Array3D.convertArrays2Arrays3D([a])
c = Converter.Array3D.convertArrays3D2Arrays(b)
T.testA(c, 1)
