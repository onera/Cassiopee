# - curve (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

# User definition of parametric curve
def f(t):
    x = t; y = t*t+1; z = 0.
    return (x,y,z)

a = D.curve(f)
t = C.newPyTree(['Base',1,a])
test.testT(t, 1)
test.writeCoverage(100)
