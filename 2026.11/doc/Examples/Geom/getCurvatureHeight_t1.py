# - getCurvatureHeight (array) -
import KCore.test as test
import Geom as D
import Converter as C
import Generator as G

# User definition of parametric curve
def f(t):
    x = t
    y = t*t+1
    z = 0.
    return (x,y,z)

# i-array
a1 = D.curve(f, 10)
b1 = D.getCurvatureHeight(a1)
test.testA([b1], 1)

# BAR-array
a2 = C.convertArray2Tetra(a1)
b2 = D.getCurvatureHeight( a2 )
test.testA([b2], 2)

# structured 2D
a3 = D.sphere((0,0,0), 1., 20)
b3 = D.getCurvatureHeight( a3 )
test.testA([b3], 3)

# QUAD
a4 = C.convertArray2Hexa(a3)
a4 = G.close(a4)
b4 = D.getCurvatureHeight( a4 )
test.testA([b4],4)

# TRI
a5 = C.convertArray2Hexa(a3)
a5 = G.close(a5)
b5 = D.getCurvatureHeight( a5 )
test.testA([b5],5)

# liste d arrays
b = D.getCurvatureHeight([a1,a2,a3])
test.testA(b, 6)
