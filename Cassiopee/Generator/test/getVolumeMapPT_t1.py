# - getVolumeMap (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as T

# Test 3D structure
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,3))
a = C.addBC2Zone(a, 'wall1','BCWall','jmin')
a = C.addBC2Zone(a, 'match1','BCMatch','imin',a,'imax',[1,2,3])
a = C.addBC2Zone(a, 'match2','BCMatch','imax',a,'imin',[1,2,3])
a = C.fillEmptyBCWith(a,'overlap','BCOverlap')
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
vol = G.getVolumeMap(a)
T.testT(vol, 1)

# Test 2D structure
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
vol = G.getVolumeMap(a)
T.testT(vol, 2)

# Test 1D structure
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (1,1,10))
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
vol = G.getVolumeMap(a)
T.testT(vol, 7)

# Test 1d non-structure bar
a = D.line((0,0,0), (1,0,0),11)
a2 = C.convertArray2Tetra(a)
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
vol = G.getVolumeMap(a)
T.testT(vol, 8)


# Test 2d non-structure hexa
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
vol = G.getVolumeMap(a)
T.testT(vol, 3)

# Test 3d non-structure hexa
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
vol = G.getVolumeMap(a)
T.testT(vol, 4)

# Test 2d non-structure tetra
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
vol = G.getVolumeMap(a)
T.testT(vol, 5)

# Test 3d non-structure tetra
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
vol = G.getVolumeMap(a)
T.testT(vol, 6)
