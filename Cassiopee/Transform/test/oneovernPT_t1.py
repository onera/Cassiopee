# - oneovern (pyTree) -
import Transform.PyTree as T
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Structure 2D
a = G.cart((0,0,0), (1,1,1), (10,10,1))
C._initVars(a, '{centers:celln}=({centers:CoordinateX}>5.)*{CoordinateX}')
C._initVars(a, '{Density}=3*{CoordinateX}*{CoordinateY}')
a = C.addBC2Zone(a, 'wall', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'overlap', 'BCOverlap', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2])
a = T.oneovern(a, (2,2,1))
test.testT(a,1)

# Structure 3D
a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._initVars(a, '{Density}=3*{CoordinateX}*{CoordinateY}')
C._initVars(a, '{centers:celln}=({centers:CoordinateX}>5.)*{CoordinateX}')
a = C.addBC2Zone(a, 'wall', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'overlap', 'BCOverlap', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2,3])
a1 = T.oneovern(a, (3,3,1))
test.testT(a1,2)

a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._initVars(a, '{centers:celln}=({centers:CoordinateX}>5.)*{CoordinateX}')
C._initVars(a,'{Density}=3*{CoordinateX}*{CoordinateY}')
a = C.addBC2Zone(a, 'wall','BCWall','imin')
a = C.addBC2Zone(a, 'overlap','BCOverlap','jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2,3])
a1 = T.oneovern(a, (2,2,1))
test.testT(a1,3)
#
a = G.cart((0,0,0), (1,1,1), (11,11,11))
C._initVars(a, '{centers:celln}=({centers:CoordinateX}>5.)*{CoordinateX}')
C._initVars(a,'{Density}=3*{CoordinateX}*{CoordinateY}')
a = C.addBC2Zone(a, 'wall','BCWall','imin')
a = C.addBC2Zone(a, 'overlap','BCOverlap','jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2,3])
a1 = T.oneovern(a, (3,3,1))
test.testT(a1,4)
#
a = G.cart((0,0,0), (1,1,1), (11,11,11))
C._initVars(a, '{centers:celln}=({centers:CoordinateX}>5.)*{CoordinateX}')
C._initVars(a,'{Density}=3*{CoordinateX}*{CoordinateY}')
a = C.addBC2Zone(a, 'wall','BCWall','imin')
a = C.addBC2Zone(a, 'overlap','BCOverlap','jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2,3])
a1 = T.oneovern(a, (2,2,1))
test.testT(a1,5)
