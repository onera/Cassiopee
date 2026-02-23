# - merge (array)
import Converter as C
import Generator as G
import Transform as T
import Geom as D
import KCore.test as test

def f(t,u):
    x = t+u
    y = t*t+1+u*u
    z = u
    return (x,y,z)

# surface grids
a = D.surface(f)
a = C.initVars(a, 'F', 1.)
b = T.splitSize(a, 100)
b1 = T.merge(b)
test.testA(b1,1)

# volume grids
distrib = G.cart((0.,0.,0.),(0.1,1,1),(11,1,1))
a = G.addNormalLayers(a, distrib)
b = T.splitSize(a, 5000)
b = T.merge(b)
test.testA(b,2)

a1 = G.cart((0,0,0),(1,1,1),(11,11,1))
a3 = G.cart((10,0,0),(1,1,1),(11,11,1))
a2 = T.rotate(a1,(0,0,0),(1,0,0),90.)
a2 = T.reorder(a2,(-1,2,3))
a = T.merge([a1,a2,a3],alphaRef=45.)
test.testA(a,3)

# mix types
a = D.surface(f)
a = C.initVars(a,'F',1.)
b = T.splitSize(a,100)
b[0] = C.convertArray2Tetra(b[0])
b[1] = C.convertArray2Tetra(b[1])
b1 = T.merge(b)
test.testA(b1,4)
