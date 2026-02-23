# - getCurvatureHeight (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import Generator.PyTree as G
import KCore.test as test

# User definition of parametric curve
def f(t):
    x = t
    y = t*t+1
    z = 0.
    return (x,y,z)

# i-array ouvert
a1 = D.curve(f, 10)
b1 = D.getCurvatureHeight(a1)
test.testT([b1], 1)

# BAR-array ouvert
a2 = C.convertArray2Tetra(a1)
b2 = D.getCurvatureHeight( a2 )
test.testT([b2], 2)

# i-array ferme
a1 = D.circle((0,0,0),1.)
b1 = D.getCurvatureHeight(a1)
test.testT([b1], 11)

# BAR-array ferme
a2 = D.circle((0,0,0),1.)
b2 = D.getCurvatureHeight(a2)
test.testT([b2], 21)

# structured 2D
a3 = D.sphere((0,0,0), 1., 20)
b3 = D.getCurvatureHeight( a3 )
test.testT([b3], 3)

# QUAD
a4 = C.convertArray2Hexa(a3)
a4 = G.close(a4)
b4 = D.getCurvatureHeight( a4 )
test.testT([b4],4)

# TRI
a5 = C.convertArray2Tetra(a3)
a5 = G.close(a5)
b5 = D.getCurvatureHeight( a5 )
test.testT([b5],5)

# liste de zones
b = D.getCurvatureHeight([a3,a4,a5])
test.testT(b, 6)

# sur un arbre
t = C.newPyTree(['BASE1',2,'BASE2',2]); t[2][1][2] = [a4]; t[2][2][2] = [a5]
t = D.getCurvatureHeight(t)
test.testT(t, 7)
