# - reorder (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import KCore.test as test
import Post.PyTree as P
import Geom.PyTree as D

# tests sur une zone de l'arbre
# reorder a Cartesian structured mesh
a = G.cart((0,0,0), (1,1,1), (8,9,20))
C._initVars(a, '{Density}={CoordinateX}**2+2*{CoordinateY}')
C._initVars(a, '{centers:F}={centers:CoordinateX}')
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imax')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmax')
a = T.reorder(a, (2,1,-3))
test.testT(a,1)

# reorder a TRI mesh
a = G.cartTetra((0,0,0), (1,1,1), (8,9,1))
C._initVars(a, '{Density}={CoordinateX}**2+2*{CoordinateY}')
C._initVars(a, '{centers:F}={centers:CoordinateX}')
a = T.reorder(a, (-1,))
test.testT(a,2)

# reorder a QUAD mesh
a = G.cartHexa((0,0,0), (1,1,1), (8,9,1))
C._initVars(a, '{Density}={CoordinateX}**2+2*{CoordinateY}')
C._initVars(a, '{centers:F}={centers:CoordinateX}')
a = T.reorder(a, (-1,))
test.testT(a,3)

# reorder a NGON mesh: sphere case
a = D.sphere6((0,0,0), 1., 20)
a = C.convertArray2NGon(a)
a = T.join(a)
a = G.close(a)
a = T.reorder(a, (1,))
G._getNormalMap(a)
a = C.initVars(a, '{centers:DensityX}=12.*{centers:sx}')
C._initVars(a, '{centers:DensityY}=12.*{centers:sy}')
C._initVars(a, '{centers:DensityZ}=12.*{centers:sz}')
resX = P.integ(a, var='centers:DensityX')
resY = P.integ(a, var='centers:DensityX')
resZ = P.integ(a, var='centers:DensityX')
test.testT(a,4)

# reorder a NGON mesh v2
a = G.cartNGon((0,0,0), (1,1,1), (8,9,1))
C._initVars(a, '{Density}={CoordinateX}**2+2*{CoordinateY}')
C._initVars(a, '{centers:F}={centers:CoordinateX}')
a = T.reorder(a, (-1,))
#test.testT(a,5)
