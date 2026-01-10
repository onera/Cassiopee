# - addkplane (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import KCore.test as test

# traitement par zone
# 3D structure + champs en noeuds + champs en centres + CL
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(11,10,2))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2,3])
a = C.addBC2Zone(a, 'nearmatch1', 'BCNearMatch', 'jmax', a, 'jmin', [1,2,3])
a = C.initVars(a, 'F', 2.); a = C.initVars(a, 'centers:G', 1.)
a = T.addkplane(a)
test.testT(a, 1)

# 2D structure + champs en noeuds + [champs en centres] + CL
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(11,10,1))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2])
a = C.addBC2Zone(a, 'nearmatch1', 'BCNearMatch', 'jmax', a, 'jmin', [1,2,3])
a = C.initVars(a, 'F', 2.); a = C.initVars(a, 'centers:G', 1.)
a = T.addkplane(a)
test.testT(a, 2)

# BAR
a = G.cart((0.,0.,0.),(0.1,1,1),(11,1,1))
a = C.convertArray2Tetra(a)
a = C.initVars(a, 'F',2.); a = C.initVars(a, 'centers:G',1.)
a = T.addkplane(a)
test.testT(a, 3)

# QUAD
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(11,10,1))
a = C.convertArray2Hexa(a)
a = C.initVars(a, 'F',2.); a = C.initVars(a, 'centers:G',1.)
a = T.addkplane(a)
test.testT(a, 4)

# TRI
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(11,10,1))
a = C.convertArray2Tetra(a)
a = C.initVars(a, 'F',2.); a = C.initVars(a, 'centers:G',1.)
a = T.addkplane(a)
test.testT(a, 5)

# 2D ME
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (5,10,1))
b = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.2), (5,10,1))
c = G.cartTetra((0.8,0.,0.), (0.1,0.1,0.2), (5,10,1))
d = G.cartTetra((0.4,-0.9,0.), (0.1,0.1,0.2), (5,10,1))
a = C.mergeConnectivity([a, b, c, d], None)
a = C.initVars(a, 'F',2.); a = C.initVars(a, 'centers:G',1.)
a = T.addkplane(a)
C.convertPyTree2File(a, 'out1.cgns')
test.testT(a, 6)
