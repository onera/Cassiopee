# - surface (PyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

# User definition of parametric surface
def f(t,u):
    x = t+u; y = t*t+1+u*u; z = u
    return (x,y,z)

a = D.surface(f)
t = C.newPyTree(['Base',2,a])
test.testT(t, 1)
test.writeCoverage(100)
