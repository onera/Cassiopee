# - curve (array) -
import Geom as D
import KCore.test as test

# User definition of parametric curve by a function
def f(t):
    x = t
    y = t*t+1
    z = 0.
    return (x,y,z)

a = D.curve(f)
test.testA([a],1)

# Definition by equation
a = D.curve('{x}=cos(2*pi*{t}); {y}=sin(2*pi*{t}); {z} = 0.')
test.testA([a],2)
test.writeCoverage(100)
