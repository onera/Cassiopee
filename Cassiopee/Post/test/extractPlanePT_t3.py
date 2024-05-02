# - extractPlane (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Geom.PyTree as D
import KCore.test as test

# test structure
m = D.sphere((0.,0.,0.),1.,N=20)
m = C.addBC2Zone(m, 'wall', 'BCWall', 'kmin')
m = C.initVars(m, 'F', 1.); m = C.initVars(m, 'centers:G', 3.)
t = C.newPyTree(['Base',2]); t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t[2][1][2].append(m)
p = P.extractPlane(t, (0.,0.,1.,-0.5))
test.testT(p, 1)

# test non structure tri
m2 = C.convertArray2Tetra(m)
m2 = C.initVars(m2, 'F', 1.); m2 = C.initVars(m2, 'centers:G', 3.)
t = C.newPyTree(['Base',2]); t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t[2][1][2].append(m2)
p = P.extractPlane(t, (0.,0.,1.,-0.5))
test.testT(p,2)

# test non structure quad
m2 = C.convertArray2Hexa(m)
m2 = C.initVars(m2, 'F', 1.); m2 = C.initVars(m2, 'centers:G', 3.)
t = C.newPyTree(['Base',2]); t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t[2][1][2].append(m2)
p = P.extractPlane(t, (0.,0.,1.,-0.5))
test.testT(p,3)
