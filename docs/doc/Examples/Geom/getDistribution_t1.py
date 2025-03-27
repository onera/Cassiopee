# - getDistribution -
import Geom as D
import KCore.test as test
import Transform as T

# User definition of parametric curve
def f(t):
    x = t
    y = t*t+1
    z = 0.
    return (x,y,z)

def g(t):
    x = 2*t
    y = t*t + 4
    z = 1.
    return (x,y,z)

# test structure
a1 = D.curve(f)
b1 = D.getDistribution(a1)
test.testA([b1],1)

# Test sur une liste
a = D.circle( (0,0,0), 1, N=100)
b = T.subzone(a, (1,1,1), (50,1,1) )
c = T.subzone(a, (50,1,1), (100,1,1) )
c2 = D.getDistribution([b,c])
test.testA(c2,3)
