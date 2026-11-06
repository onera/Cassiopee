# - initVars (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Create a function
def F(x1, x2): return 3.*x1+2.*x2

# Sur une zone
a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.initVars(a, 'centers:celln', 2.)
a = C.initVars(a, 'Density', F, ['CoordinateX','CoordinateY'])
a = C.initVars(a, 'centers:G', F, ['centers:CoordinateX','centers:CoordinateY'])
t = C.newPyTree(['Base', a])
test.testT(t, 1)

# Avec une formule
a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.initVars(a, '{F} = 3*{CoordinateX}')
a = C.initVars(a, '{centers:G} = 4 * {centers:CoordinateX}')
t = C.newPyTree(['Base', a])
test.testT(t, 11)

# Sur une zone NGon
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
a = C.initVars(a, 'centers:celln', 2.)
a = C.initVars(a,'{Density}=3*{CoordinateX}+2*{CoordinateY}')
#a = C.initVars(a,'centers:G', F,['centers:CoordinateX','centers:CoordinateY'])
t = C.newPyTree(['Base', a])
test.testT(t, 12)

# Sur un arbre
a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((10,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base', a,b])
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.initVars(t, 'centers:celln', 2.)
t = C.initVars(t, 'Density', F, ['CoordinateX','CoordinateY'])
test.testT(t, 2)

# Sur une liste de zones
A = [a,b]
A = C.initVars(A, 'centers:celln', 2.)
A = C.initVars(A, 'Density', F, ['CoordinateX','CoordinateY'])
t = C.newPyTree(['Base']); t[2][1][2] += A
test.testT(t, 3)
