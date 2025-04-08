# - rmBCDataVars (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Sur une zone
a = G.cart((0,0,0),(1,1,1),(10,10,10))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'wall2', 'BCWall', 'jmax')
C._initBCDataSet(a,'{var1}=1.')
C._initBCDataSet(a,'{var2}=2.')
C._initBCDataSet(a,'{var3}=3.')
a = C.rmBCDataVars(a,'var1')
a = C.rmBCDataVars(a,['var2','var3'])

t = C.newPyTree(['Base',a])
test.testT(t, 1)

# Sur un arbre
a = G.cart((0,0,0),(1,1,1),(10,10,10))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'wall2', 'BCWall', 'jmax')

b = G.cart((10,0,0),(1,1,1),(10,10,10))
b = C.addBC2Zone(b, 'bcsym1', 'BCSymmetryPlane', 'imin')
b = C.addBC2Zone(b, 'bcsym2', 'BCSymmetryPlane', 'imax')

t = C.newPyTree(['Base',a,b])

C._initBCDataSet(t,'{var1}=1.')
C._initBCDataSet(t,'{var2}=2.')
C._initBCDataSet(t,'{var3}=3.')

t = C.rmBCDataVars(t,'var1')
t = C.rmBCDataVars(t,['var2','var3'])

test.testT(t, 2)

# Sur une liste de zones
A = C.rmBCDataVars([a,b], 'var1')
t = C.newPyTree(['Base']); t[2][1][2] += A
test.testT(t, 3)

B = C.rmBCDataVars([a,b], ['var2','var3'])
t = C.newPyTree(['Base']); t[2][1][2] += B
test.testT(t, 4)
