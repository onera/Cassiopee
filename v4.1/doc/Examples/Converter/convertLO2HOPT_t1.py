# - convertLO2HO (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

# BAR -> BAR_3
a = D.line((0,0,0), (1,0,0), N=3)
a = C.convertArray2Hexa(a)
a = C.convertLO2HO(a, 0)
test.testT(a, 1)

# TRI -> TRI_6
a = D.triangle((0,0,0), (1,0,0), (1,1,0))
a = C.convertLO2HO(a, 0)
test.testT(a, 2)

# QUAD -> QUAD_8
a = D.quadrangle((0,0,0), (1,0,0), (1,1,0), (0,1,0))
a = C.convertLO2HO(a, 0)
test.testT(a, 3)

# QUAD -> QUAD_9
a = D.quadrangle((0,0,0), (1,0,0), (1,1,0), (0,1,0))
a = C.convertLO2HO(a, 1)
test.testT(a, 4)

# TETRA -> TETRA_10
a = G.cart((0,0,0), (1,1,1), (2,2,2))
a = C.convertArray2Tetra(a)
a = T.subzone(a, [0], type='elements')
a = C.convertLO2HO(a, 0)
test.testT(a, 5)
