# - surface (array) -
import Geom as D
import KCore.test as test

# User definition of parametric surface by function
def f(t,u):
    x = t+u
    y = t*t+1+u*u
    z = u
    return (x,y,z)

a = D.surface(f)
test.testA([a],1)

# Definition by formula
a = D.surface('{x} = cos(pi*{t}); {y} = sin(pi*{u}); {z} = {t}*{u}')
test.testA([a],2)

test.writeCoverage(100)
