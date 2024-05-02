# - join (pyTree) -
import Generator.PyTree as G
import Geom.PyTree as D
import Transform.PyTree as T
import Converter.PyTree as C
import KCore.test as test

# Join 2 NON-STRUCT TETRA
a1 = G.cartTetra((0.,0.,0.), (1.,1.,1), (11,11,10))
a2 = G.cartTetra((10.,0.,0.), (1.,1.,1), (10,10,10))
C._initVars(a1, 'F', 2); C._initVars(a1, 'centers:G', 1); C._initVars(a1, 'centers:G2', 10)
C._initVars(a2, 'F', 3); C._initVars(a2, 'centers:G2', 10); C._initVars(a2, 'centers:G', 3)
t = C.newPyTree(['Base',3])
a = T.join(a1, a2); t[2][1][2].append(a)
test.testT(t, 1)

# mix struct 3D + HEXA
a1 = G.cart((0.,0.,0.), (1.,1.,1), (11,11,10))
C._addBC2Zone(a1,'wall','BCWall','imin')
C._addBC2Zone(a1,'ov','BCOverlap','imax')
C._addBC2Zone(a1,'match1','BCMatch','jmin',a1, 'jmax')
a2 = G.cartHexa((10.,0.,0.), (1.,1.,1), (11,11,10))
C._initVars(a1, 'F', 2); C._initVars(a1, 'centers:G', 1)
C._initVars(a2, 'F', 3); C._initVars(a2, 'centers:G', 3)
t = C.newPyTree(['Base',3])
a = T.join(a1, a2); t[2][1][2].append(a)
test.testT(t, 2)

# Join 2 NON-STRUCT PENTA
a1 = G.cartPenta((0.,0.,0.), (1.,1.,1), (11,11,10))
a2 = G.cartPenta((10.,0.,0.), (1.,1.,1), (10,10,10))
a1 = C.initVars(a1, 'F', 2); a1 = C.initVars(a1, 'centers:G', 1)
a2 = C.initVars(a2, 'F', 3); a2 = C.initVars(a2, 'centers:G', 3)
t = C.newPyTree(['Base',2])
a = T.join(a1, a2); t[2][1][2].append(a)
test.testT(t, 9)

# Join 2 NON-STRUCT PYRA
a1 = G.cartPyra((0.,0.,0.), (1.,1.,1), (11,11,10))
a2 = G.cartPyra((10.,0.,0.), (1.,1.,1), (10,10,10))
a1 = C.initVars(a1, 'F', 2); a1 = C.initVars(a1, 'centers:G', 1)
a2 = C.initVars(a2, 'F', 3); a2 = C.initVars(a2, 'centers:G', 3)
t = C.newPyTree(['Base',2])
a = T.join(a1, a2); t[2][1][2].append(a)
test.testT(t, 10)

# Join 2 NON-STRUCT QUADS
a1 = G.cartTetra((0.,0.,0.), (1.,1.,1), (11,11,1))
a2 = G.cartTetra((10.,0.,0.), (1.,1.,1), (10,10,1))
a1 = C.initVars(a1, 'F', 2); a1 = C.initVars(a1, 'centers:G', 1)
a2 = C.initVars(a2, 'F', 3); a2 = C.initVars(a2, 'centers:G', 3)
t = C.newPyTree(['Base',2])
a = T.join(a1, a2); t[2][1][2].append(a)
test.testT(t, 3)

# Join 1 STRUCT et 1 NON-STRUCT QUAD
a1 = G.cart((13.,0.,0.), (1.,1.,1), (11,11,1))
a1 = C.addBC2Zone(a1,'wall','BCWall','imin')
a1 = C.addBC2Zone(a1,'ov','BCOverlap','imax')
a1 = C.addBC2Zone(a1,'match1','BCMatch','jmin',a1, 'jmax')
a2 = G.cartHexa((9.,0.,0.), (1.,1.,1), (10,10,1))
a1 = C.initVars(a1, 'F', 2); a1 = C.initVars(a1, 'centers:G', 1)
a2 = C.initVars(a2, 'F', 3); a2 = C.initVars(a2, 'centers:G', 3)
t = C.newPyTree(['Base',2])
a = T.join(a1, a2); t[2][1][2].append(a)
test.testT(t, 4)

# Join 2 STRUCT-1D
a1 = D.line((0.,0.,0.), (1.,0.,0), 100)
a2 = D.line((1.,0.,0.), (1.,1,0), 100)
a1 = C.addBC2Zone(a1,'wall','BCWall','imin')
a1 = C.addBC2Zone(a1,'ov','BCOverlap','imax')
a1 = C.initVars(a1, 'F', 2); a1 = C.initVars(a1, 'centers:G', 1)
a2 = C.initVars(a2, 'F', 3); a2 = C.initVars(a2, 'centers:G', 3)
t = C.newPyTree(['Base',2])
a = T.join(a1, a2); t[2][1][2].append(a)
test.testT(t, 5)

# Join 2 BARS
a1 = D.line((0.,0.,0.), (1.,0.,0), 100)
a2 = D.line((1.,0.,0.), (1.,1,0), 100)
a1 = C.convertArray2Tetra(a1)
a2 = C.convertArray2Tetra(a2)
a1 = C.initVars(a1, 'F', 2); a1 = C.initVars(a1, 'centers:G', 1)
a2 = C.initVars(a2, 'F', 3); a2 = C.initVars(a2, 'centers:G', 3)
t = C.newPyTree(['Base',1])
a = T.join(a1, a2); t[2][1][2].append(a)
test.testT(t, 6)

# Join 1 BAR et 1 i-array
a1 = D.line((0.,0.,0.), (1.,0.,0), 100)
a1 = C.addBC2Zone(a1,'wall','BCWall','imin')
a1 = C.addBC2Zone(a1,'ov','BCOverlap','imax')
a2 = D.line((1.,0.,0.), (1.,1,0), 100)
a2 = C.convertArray2Tetra(a2)
a1 = C.initVars(a1, 'F', 2); a1 = C.initVars(a1, 'centers:G', 1)
a2 = C.initVars(a2, 'F', 3); a2 = C.initVars(a2, 'centers:G', 3)
t = C.newPyTree(['Base',1])
a = T.join(a1, a2); t[2][1][2].append(a)
test.testT(t, 7)

# Join 2  STRUCT 3D
a1 = G.cart((0.,0.,0.), (1.,1.,1), (11,11,10))
a2 = G.cart((10.,0.,0.), (1.,1.,1), (11,11,10))
a1 = C.addBC2Zone(a1,'match1','BCMatch','imax',a2,'imin',[1,2,3])
a2 = C.addBC2Zone(a2,'match1','BCMatch','imin',a1,'imax',[1,2,3])
a1 = C.initVars(a1, 'F', 2); a1 = C.initVars(a1, 'centers:G', 1)
a2 = C.initVars(a2, 'F', 3); a2 = C.initVars(a2, 'centers:G', 3)
t = C.newPyTree(['Base']); a = T.join(a1, a2); t[2][1][2].append(a)
test.testT(t,8)
