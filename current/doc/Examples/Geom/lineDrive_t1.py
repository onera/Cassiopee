# - lineDrive (array) -
import Geom as D
import Generator as G
import KCore.test as test
import Converter as C

# 1D structure
a = D.naca(12.)
b = D.line((0,0,0), (0,0.,1.))
c = D.lineDrive(a, b)
test.testA([c], 1)

# 1D structure + champ
a = D.circle((0,0,0), 1)
a = C.addVars(a, 'var')
b = D.line((0,0,0), (0,0.,1.))
c = D.lineDrive(a, b)
test.testA([c], 2)

# 2D structure
a = G.cylinder((0,0,0), 1, 2, 360, 0, 1, (50,21,1))
a = C.addVars(a, 'var')
b = D.line((0,0,0), (0,0.,1.))
c = D.lineDrive(a, b)
test.testA([c], 3)

# BAR
a = D.line((0,0,0), (1,0,0), N=10)
a = C.convertArray2Tetra(a)
b = D.line((0,0,0), (0,0.,1.), N=10)
c = D.lineDrive(a, b)
test.testA([c], 4)

# QUAD
a = G.cart((0,0,0), (1,1,1), (10,10,1))
a = C.convertArray2Hexa(a)
b = D.line((0,0,0), (0,0.,1.), N=10)
c = D.lineDrive(a, b)
test.testA([c], 5)

# TRI
a = G.cart((0,0,0), (1,1,1), (10,10,1))
a = C.convertArray2Tetra(a)
b = D.line((0,0,0), (0,0.,1.), N=10)
c = D.lineDrive(a, b)
test.testA([c], 6)
