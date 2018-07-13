# - cpVars (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0),(1,1,1),(10,10,10))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
b = G.cart((0,0,0),(1,1,1),(10,10,10)); b[0] = 'cart2'
c = G.cart((0,0,0),(1,1,1),(10,10,10)); c[0] = 'cart3'
c = C.initVars(c,'F',10.); c = C.initVars(c,'centers:G',20.)

# Copie variables a -> b (nodes) avec creation
b = C.cpVars(a, 'F', b, 'F')
# Copie variables a -> c (nodes) sans creation
c = C.cpVars(a, 'F', c, 'F')
t = C.newPyTree(['Base']); t[2][1][2] += [b,c]
test.testT(t, 1)

# Copie variables a -> b (centers) avec creation
b = C.cpVars(a, 'centers:G', b, 'centers:G')
# Copie variables a -> c (centers) sans creation
c = C.cpVars(a, 'centers:G', c, 'centers:G')
t = C.newPyTree(['Base']); t[2][1][2] += [b,c]
test.testT(t, 2)

a = G.cart((0,0,0),(1,1,1),(10,10,10))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
b = G.cart((0,0,0),(1,1,1),(10,10,10)); b[0] = 'cart2'
c = G.cart((0,0,0),(1,1,1),(10,10,10)); c[0] = 'cart3'
c = C.initVars(c,'F',10.); c = C.initVars(c,'centers:G',20.)
b = C.node2Center(b)
c = C.node2Center(c)

# Copie variables a -> b (nodes/centers) avec creation
b = C.cpVars(a, 'centers:G', b, 'G')
# Copie variables a -> c (nodes/centers) sans creation
c = C.cpVars(a, 'centers:G', c, 'G')

t = C.newPyTree(['Base']); t[2][1][2] += [b,c]
test.testT(t, 3)
