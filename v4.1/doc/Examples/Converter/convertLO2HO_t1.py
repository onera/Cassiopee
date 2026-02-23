# - convertLO2HO (array) -
import Converter as C
import Geom as D
import Generator as G
import Transform as T
import KCore.test as test

############# Order 2 ################

# BAR -> BAR_3
a = D.line((0,0,0), (1,0,0), N=3)
a = C.convertArray2Hexa(a)
a = C.convertLO2HO(a, 0, 2)
test.testA(a, 1)

# TRI -> TRI_6
a = D.triangle((0,0,0), (1,0,0), (1,1,0))
a = C.convertLO2HO(a, 0, 2)
test.testA(a, 2)

# QUAD -> QUAD_8
a = D.quadrangle((0,0,0), (1,0,0), (1,1,0), (0,1,0))
a = C.convertLO2HO(a, 0, 2)
test.testA(a, 3)

# QUAD -> QUAD_9
a = D.quadrangle((0,0,0), (1,0,0), (1,1,0), (0,1,0))
a = C.convertLO2HO(a, 1, 2)
test.testA(a, 4)

# TETRA -> TETRA_10
a = G.cart((0,0,0), (1,1,1), (2,2,2))
a = C.convertArray2Tetra(a)
a = T.subzone(a, [0], type='elements')
a = C.convertLO2HO(a, 0, 2)
test.testA(a, 5)

############# Order 3 ################

# BAR -> BAR_4
a = D.line((0,0,0), (1,0,0), N=3)
a = C.convertArray2Hexa(a)
a = C.convertLO2HO(a, 0, 3)
test.testA(a, 6)

# TRI -> TRI_9
a = D.triangle((0,0,0), (1,0,0), (1,1,0))
a = C.convertLO2HO(a, 0, 3)
test.testA(a, 7)

# TRI -> TRI_10
a = D.triangle((0,0,0), (1,0,0), (1,1,0))
a = C.convertLO2HO(a, 1, 3)
test.testA(a, 8)

# QUAD -> QUAD_12
a = D.quadrangle((0,0,0), (1,0,0), (1,1,0), (0,1,0))
a = C.convertLO2HO(a, 0, 3)
test.testA(a, 9)

# QUAD -> QUAD_16
a = D.quadrangle((0,0,0), (1,0,0), (1,1,0), (0,1,0))
a = C.convertLO2HO(a, 1, 3)
test.testA(a, 10)

############# Order 4 ################

# TRI -> TRI_12
a = D.triangle((0,0,0), (1,0,0), (1,1,0))
a = C.convertLO2HO(a, 0, 4)
test.testA(a, 11)

# TRI -> TRI_15
a = D.triangle((0,0,0), (1,0,0), (1,1,0))
a = C.convertLO2HO(a, 1, 4)
test.testA(a, 12)

# QUAD -> QUAD_P4_16 ######## SE PLANTE DANS GETFROMARRAY ????
#a = D.quadrangle((0,0,0), (1,0,0), (1,1,0), (0,1,0))
#a = C.convertLO2HO(a, 0, 4)
#test.testA(a, 13)

# QUAD -> QUAD_25
a = D.quadrangle((0,0,0), (1,0,0), (1,1,0), (0,1,0))
a = C.convertLO2HO(a, 1, 4)
test.testA(a, 14)
