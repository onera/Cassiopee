# - getCurvatureAngle (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

a = D.polyline([(0.,0.,0.),(1.,1.,0.),(2.,0.,0.)])
a = D.getCurvatureAngle( a )
t = C.newPyTree(['Base',1]); t[2][1][2].append(a)
test.testT(t, 1)

# User definition of parametric curve
def f(t):
    x = t; y = t*t+1; z = 0.
    return (x,y,z)

# i-array
a1 = D.curve(f, 10)
b1 = D.getCurvatureAngle(a1)
test.testT([b1], 2)

# BAR-array
a2 = C.convertArray2Tetra(a1)
b2 = D.getCurvatureAngle(a2)
t = C.newPyTree(['Base',1]); t[2][1][2].append(b2)
test.testT(t, 3)

# TRI-array
#a3 = D.sphere((0,0,0), 1., 20)
#a3 = C.convertArray2Tetra(a3)
#a3 = G.close(a3)
#b3 = D.getCurvatureAngle( a3 )
#t = C.newPyTree(['Base',3]); t[2][1][2].append(b3)
#test.testT(t, 4)
