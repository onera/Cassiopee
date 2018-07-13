# - getNormalMap (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# Test 2D structure
a = D.sphere((0,0,0), 1, 50)
a = C.initVars(a,'Density',2.); a = C.initVars(a,'centers:cellN',1.)
n = G.getNormalMap(a)
test.testT(n,1)

# Test 2d non-structure quad
n = C.convertArray2Hexa(a)
n = C.initVars(n,'Density',2.); n = C.initVars(n,'centers:cellN',1.)
n = G.getNormalMap(n)
test.testT(n,2)

# Test 2d non-structure tri
n = C.convertArray2Tetra(a)
n = C.initVars(n,'Density',2.); n = C.initVars(n,'centers:cellN',1.)
n = G.getNormalMap(n)
test.testT(n,3)


# test NGon surfacique
a = D.sphere((0,0,0), 1, 10)
a = C.convertArray2NGon(a)
a = C.initVars(a,'Density',2.); a = C.initVars(a,'centers:cellN',1.)
a = G.getNormalMap(a)
test.testT(a,4)
