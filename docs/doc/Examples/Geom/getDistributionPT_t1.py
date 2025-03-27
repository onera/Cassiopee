# - getDistribution (pyTree) -
import Geom.PyTree as D
import KCore.test as test

a = D.line((0.,0.,0.), (1.,0.,0), 100)
a = D.getDistribution(a)
test.testT(a, 1)

# User definition of parametric curve
def f(t):
    x = t; y = t*t+1; z = 0.
    return (x,y,z)

def g(t):
    x = 2*t; y = t*t + 4; z = 1.
    return (x,y,z)

# test structure
a1 = D.curve(f)
b1 = D.getDistribution(a1)
test.testT(b1, 2)
