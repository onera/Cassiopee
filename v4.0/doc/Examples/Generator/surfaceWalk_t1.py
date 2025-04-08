# - surfaceWalk (array)
import Geom as D
import Transform as T
import Generator as G
import KCore.test as test

# User definition of parametric curve
def f(t,u):
    x = t+u; y = t*t+1+u*u; z = u
    return (x,y,z)

# Array definition of geometry
a = D.surface(f)
c = D.circle((1.2,1.7,0.6), 0.1)
c = T.rotate(c, (1.2,1.7,0.6), (0,1,0), 90.)
c = T.reorder(c,(-1,2,3))
c = T.projectOrtho(c,[a])
h = G.cart((0.,0.,0.),(0.01,1,1),(10,1,1))
r = G.surfaceWalk([a], c, h, niter=100)
test.testA([r])
