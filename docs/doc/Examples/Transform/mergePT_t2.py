# - merge (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Geom.PyTree as D
import KCore.test as test
def f(t,u):
    x = t+u
    y = t*t+1+u*u
    z = u
    return (x,y,z)

# surface grids
a = D.surface(f)
t = C.newPyTree(['Base', 2])
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
t[2][1][2] += b
t = C.initVars(t, 'F', 1.); t = C.initVars(t, 'centers:G', 2.)
res = T.merge(t, dir=2)
t[2][1][2] = res; test.testT(t,1)
res = T.merge(t,dir=1)
t[2][1][2] = res; test.testT(t,2)

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
        b.append(T.subzone(a,(i1,j1,1),(i2,j2,10)))
        j1 = j2
    i1 = i2
t = C.newPyTree(['Base']); t[2][1][2]+=b
t = C.addBC2Zone(t,'wall','BCWall','kmin')
t = C.addBC2Zone(t,'overlap','BCOverlap','kmax')
t = C.initVars(t,'F',1.); t = C.initVars(t,'centers:G',2.)
res = T.merge(t, dir=2)
t[2][1][2] = res; test.testT(t,3)
res = T.merge(t, dir=1)
t[2][1][2] = res; test.testT(t,4)
res = T.merge(t, dir=3)
t[2][1][2] = res; test.testT(t,5)

a1 = G.cart((0,0,0),(1,1,1),(11,11,1))
a3 = G.cart((10,0,0),(1,1,1),(11,11,1))
a2 = T.rotate(a1,(0,0,0),(1,0,0),90.)
a2 = T.reorder(a2,(-1,2,3))
t = C.newPyTree(['Base', 2]); t[2][1][2] += [a1,a2,a3]
t = C.addBC2Zone(t,'wall','BCWall','imin')
t = C.addBC2Zone(t,'overlap','BCOverlap','imax')
t = C.initVars(t,'F',1.); t = C.initVars(t,'centers:G',2.)
res = T.merge(t, dir=1, alphaRef=45.)
t[2][1][2] = res; test.testT(t,6)

t = C.newPyTree(['Base',2]); t[2][1][2]+=[a1,a2,a3]
t = C.addBC2Zone(t,'wall','BCWall','kmin')
t = C.addBC2Zone(t,'overlap','BCOverlap','kmax')
t = C.initVars(t,'F',1.); t = C.initVars(t,'centers:G',2.)
res = T.merge(t,dir=2,alphaRef=45.)
t[2][1][2] = res; test.testT(t,7)

t = C.newPyTree(['Base',2]); t[2][1][2]+=[a1,a2,a3]
t = C.addBC2Zone(t,'wall','BCWall','kmin')
t = C.addBC2Zone(t,'overlap','BCOverlap','kmax')
t = C.initVars(t,'F',1.); t = C.initVars(t,'centers:G',2.)
res = T.merge(t,dir=2)
t[2][1][2] = res; test.testT(t,8)
