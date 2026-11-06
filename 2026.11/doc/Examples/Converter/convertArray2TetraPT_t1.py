# - convertArray2Tetra (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Sur une zone structuree 3D + champ noeuds + champ en centres
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,5))
a = C.initVars(a, '{F}={CoordinateX}')
a = C.initVars(a, '{centers:G}={centers:CoordinateY}')
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.convertArray2Tetra(a)
t = C.newPyTree(['Base',a])
test.testT(t, 1)

# Sur un arbre
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,5))
b = G.cart((10,0.,0.), (0.1,0.1,0.2), (10,10,5))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
t = C.newPyTree(['Base',a,b])
t = C.convertArray2Tetra(t)
test.testT(t, 2)

# Sur une liste de zones
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,5))
b = G.cart((10,0.,0.), (0.1,0.1,0.2), (10,10,5))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
A = C.convertArray2Tetra([a,b])
t = C.newPyTree(['Base']); t[2][1][2] = t[2][1][2] + A
test.testT(t, 3)

# Sur un HEXA
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,5))
a = C.initVars(a, '{F}={CoordinateX}')
a = C.initVars(a, '{centers:G}={centers:CoordinateY}')
a = C.convertArray2Tetra(a)
t = C.newPyTree(['Base',a])
test.testT(t, 4)

# Sur un QUAD
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a = C.initVars(a, '{F}={CoordinateX}')
a = C.initVars(a, '{centers:G}={centers:CoordinateY}')
a = C.convertArray2Tetra(a)
t = C.newPyTree(['Base',2]); t[2][1][2].append(a)
test.testT(t, 5)
