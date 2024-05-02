# - merge (array) -
import Converter as C
import Generator as G
import Transform as T
import Geom as D
import KCore.test as test

# tests with directions
def f(t,u):
    x = t+u
    y = t*t+1+u*u
    z = u
    return (x,y,z)
# surface grids
a = D.surface(f)
a = C.initVars(a,'F',1.)
b = []
i1 = 1
for i in range(10):
    i2 = (i+1)*10
    j1 = 1
    for j in range(10):
        j2 = (j+1)*10
        b.append(T.subzone(a,(i1,j1,1),(i2,j2,1)))
        j1 = j2
    i1 = i2
b1 = T.merge(b,dir=2)
test.testA(b1,1)
b1 = T.merge(b1,dir=1)
test.testA(b1,2)

# volume grids
distrib = G.cart((0.,0.,0.),(0.1,1,1),(10,1,1))
a = G.addNormalLayers(a,distrib)
i1 = 1
b = []
for i in range(10):
    i2 = (i+1)*10
    j1 = 1
    for j in range(10):
        j2 = (j+1)*10
        k1 = 1
        for k in range(2):
            k2 = (k+1)*5
            b.append(T.subzone(a,(i1,j1,k1),(i2,j2,k2)))
            k1 = k2
        j1 = j2
    i1 = i2
b1 = T.merge(b,dir=2)
test.testA(b1,3)
b1 = T.merge(b1,dir=1)
test.testA(b1,4)
b1 = T.merge(b,dir=3)
test.testA(b1,5)
#
a1 = G.cart((0,0,0),(1,1,1),(11,11,1))
a3 = G.cart((10,0,0),(1,1,1),(11,11,1))
a2 = T.rotate(a1,(0,0,0),(1,0,0),90.)
a2 = T.reorder(a2,(-1,2,3))
a = T.merge([a1,a2,a3],dir=1,alphaRef=45.)
test.testA(a,6)
a = T.merge([a1,a2,a3],dir=2,alphaRef=45.)
test.testA(a,7)
a = T.merge([a1,a2,a3],dir=2)
test.testA(a,8)
