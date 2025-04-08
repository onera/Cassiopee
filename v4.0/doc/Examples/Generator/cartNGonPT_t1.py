# - cartNGon (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# 1D
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (11,1,1))
t = C.newPyTree(['Base']); t[2][1][2].append(a)
test.testT(t,1)
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (1,11,1))
t = C.newPyTree(['Base']); t[2][1][2].append(a)
test.testT(t,2)
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (1,1,11))
t = C.newPyTree(['Base']); t[2][1][2].append(a)
test.testT(t,3)
# 2D
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (11,11,1))
t = C.newPyTree(['Base']); t[2][1][2].append(a)
test.testT(t,4)
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (11,1,11))
t = C.newPyTree(['Base']); t[2][1][2].append(a)
test.testT(t,5)
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (1,11,11))
t = C.newPyTree(['Base']); t[2][1][2].append(a)
test.testT(t,6)
# 3D
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (11,11,11))
t = C.newPyTree(['Base']); t[2][1][2].append(a)
test.testT(t,7)
test.writeCoverage(100)
