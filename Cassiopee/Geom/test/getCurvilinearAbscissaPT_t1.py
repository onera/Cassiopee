# - getCurvilinearAbscissa (pyTree)-
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

a = D.line((0.,0.,0.), (1.,0.,0), 100)
a = D.getCurvilinearAbscissa( a )
t = C.newPyTree(['Base',1]); t[2][1][2].append(a)
test.testT(t, 1)

# User definition of parametric curve
def f(t):
    x = t; y = t*t+1; z = 0.
    return (x,y,z)

def g(t):
    x = 2*t; y = t*t + 4; z = 1.
    return (x,y,z)

# test structure
a1 = D.curve(f)
b1 = D.getCurvilinearAbscissa( a1 )
t = C.newPyTree(['Base',1]); t[2][1][2].append(b1)
test.testT(t, 2)

# test non structure BAR
a2 = C.convertArray2Tetra(a1)
b2 = D.getCurvilinearAbscissa( a2 )
t = C.newPyTree(['Base',1]); t[2][1][2].append(b2)
test.testT(t, 3)
