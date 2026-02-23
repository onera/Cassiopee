# - rmGhostCells (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Converter.Internal as Internal
import KCore.test as test

# 3D sur une zone
a = D.sphere((0,0,0),1,N=10)
dhk = G.cart((0,0,0),(0.1,1,1),(2,1,1))
a = G.addNormalLayers(a,dhk)
a = C.convertArray2NGon(a)
a = C.initVars(a, '{centers:F}={centers:CoordinateX}')
a = C.initVars(a,'{Density}={CoordinateX}')
t = C.newPyTree(['Base']); t[2][1][2] = [a]
t = Internal.addGhostCells(t, t, 2, adaptBCs=0,fillCorner=1)
t = Internal.rmGhostCells(t, t, 2)
test.testT(t,1)

# 3D sur un arbre
a = D.sphere((0,0,0),1,N=10)
dhk = G.cart((0,0,0),(0.1,1,1),(2,1,1))
a = G.addNormalLayers(a,dhk)
a = C.convertArray2NGon(a)
a = C.initVars(a, '{centers:F}={centers:CoordinateX}**2')
a = C.initVars(a,'{Density}=3*{CoordinateX}+2*{CoordinateY}')
t = C.newPyTree(['Base']); t[2][1][2] = [a]
t = Internal.addGhostCells(t, t, 1, adaptBCs=0,fillCorner=1)
t = Internal.rmGhostCells(t, t, 1)
test.testT(t,2)

#  2D sans frontieres exterieures
a = D.sphere((0,0,0),1,N=10)
a = C.convertArray2NGon(a)
a = C.initVars(a, '{centers:F}={centers:CoordinateX}**2')
a = C.initVars(a,'{Density}=3*{CoordinateX}+2*{CoordinateY}')
t = C.newPyTree(['Base']); t[2][1][2] = [a]
t = Internal.addGhostCells(t, t, 1, adaptBCs=0,fillCorner=1)
t = Internal.rmGhostCells(t, t, 1)
test.testT(t,3)

# 2D avec frontieres exterieures
a = G.cartNGon((0,0,0),(1,1,1),(21,21,1))
a = C.initVars(a, '{centers:F}={centers:CoordinateX}**2')
a = C.initVars(a,'{Density}=3*{CoordinateX}+2*{CoordinateY}')
t = C.newPyTree(['Base']); t[2][1][2] = [a]
t = Internal.addGhostCells(t, t, 1, adaptBCs=0,fillCorner=1)
t = Internal.rmGhostCells(t, t, 1)
test.testT(t,4)
