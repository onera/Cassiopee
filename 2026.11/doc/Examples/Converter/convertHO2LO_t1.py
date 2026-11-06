# - convertHO2LO (array) -
import Converter as C
import Geom as D
import Generator as G
import Transform as T
import KCore.test as test

# BAR -> BAR_3
a = D.line((0,0,0), (1,0,0), N=3)
a = C.convertArray2Hexa(a)
a = C.convertLO2HO(a, 0)
a = C.convertHO2LO(a, 0)
test.testA(a, 1)

# TRI -> TRI_6
a = D.triangle((0,0,0), (1,0,0), (1,1,0))
a = C.convertLO2HO(a, 0)
a = C.convertHO2LO(a, 0)
test.testA(a, 2)

# QUAD -> QUAD_8
a = D.quadrangle((0,0,0), (1,0,0), (1,1,0), (0,1,0))
a = C.convertLO2HO(a, 0)
a = C.convertHO2LO(a, 0)
test.testA(a, 3)

# QUAD -> QUAD_9
a = D.quadrangle((0,0,0), (1,0,0), (1,1,0), (0,1,0))
a = C.convertLO2HO(a, 1)
a = C.convertHO2LO(a, 0)
test.testA(a, 4)

# TETRA -> TETRA_10
a = G.cart((0,0,0), (1,1,1), (2,2,2))
a = C.convertArray2Tetra(a)
a = T.subzone(a, [0], type='elements')
a = C.convertLO2HO(a, 0)
a = C.convertHO2LO(a, 0)
test.testA(a, 5)
