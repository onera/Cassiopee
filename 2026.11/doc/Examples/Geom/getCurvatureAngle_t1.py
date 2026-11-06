# - getCurvatureAngle (array) -
import KCore.test as test
import Geom as D
import Converter as C

# User definition of parametric curve
def f(t):
    x = t
    y = t*t+1
    z = 0.
    return (x,y,z)

# i-array
a1 = D.curve(f, 10)
b1 = D.getCurvatureAngle( a1 )
test.testA([b1], 1)

# BAR-array
a2 = C.convertArray2Tetra(a1)
b2 = D.getCurvatureAngle( a2 )
test.testA([b2], 2)

# TRI-array
#a3 = D.sphere((0,0,0), 1., 20)
#a3 = C.convertArray2Tetra(a3)
#a3 = G.close(a3)
#b3 = D.getCurvatureAngle( a3 )
#test.testA([b3], 3)
