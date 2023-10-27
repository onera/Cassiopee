# - zip (array) -
import Converter as C
import Generator as G
import KCore.test as test

a1 = G.cart((0,0,0), (1,1,1), (10,10,1))
a2 = G.cart((9+1.e-2,0,0), (1,1,1), (10,10,1))
a3 = G.cart((0,-5.01,0),(1,1,1),(19,6,1))
a4 = G.cart((0,9.0001,0),(1,1,1),(10,6,1))
a5 = G.cart((9.01,9.0002,0),(1,1,1),(10,6,1))

# Merge vertices at borders if they are distant from eps
zones = [a3,a1,a2,a4,a5]; zones = C.addVars(zones,'F')
B = C.copy(zones)
B = G.zip(B, 1e-1)
test.testA(B, 1)

B = C.convertArray2Hexa(zones)
B = G.zip(B, 1e-1)
test.testA(B, 2)

B = C.convertArray2Tetra(zones)
B = G.zip(B, 1e-1)
test.testA(B, 3)

B = C.convertArray2NGon(zones)
B = G.zip(B, 1e-1)
test.testA(B, 3)
