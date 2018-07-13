# - close (array) -
import Converter as C
import Generator as G
import Transform as T
import KCore.test as test

a1 = G.cart((0,0,0), (1,1,1), (10,10,1))
a2 = G.cart((9+1.e-2,0,0), (1,1,1), (10,10,1))
a3 = G.cart((0,-5.01,0),(1,1,1),(19,6,1))
a4 = G.cart((0,9.0001,0),(1,1,1),(10,6,1))
a5 = G.cart((9.01,9.0002,0),(1,1,1),(10,6,1))

# Close une liste de maillages structures
zones = [a3,a1,a2,a4,a5]; zones = C.addVars(zones,'F')
B = G.close(zones, 1e-1)
test.testA(B, 1)

