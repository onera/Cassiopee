# - getVolumeMap (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

# --- Structured grids ---
# 3D structured
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,3))
a = C.addBC2Zone(a, 'wall1','BCWall','jmin')
a = C.addBC2Zone(a, 'match1','BCMatch','imin',a,'imax',[1,2,3])
a = C.addBC2Zone(a, 'match2','BCMatch','imax',a,'imin',[1,2,3])
a = C.fillEmptyBCWith(a,'overlap','BCOverlap')
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
vol = G.getVolumeMap(a)
test.testT(vol, 1)

# 2D structured
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
vol = G.getVolumeMap(a)
test.testT(vol, 2)

# 1D structured
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (1,1,10))
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
vol = G.getVolumeMap(a)
test.testT(vol, 7)

# --- BE ---
# 1d unstructured bar
a = D.line((0,0,0), (1,0,0),11)
a2 = C.convertArray2Tetra(a)
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
vol = G.getVolumeMap(a)
test.testT(vol, 8)

# 2D unstructured quad
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
vol = G.getVolumeMap(a)
test.testT(vol, 3)

# 3D unstructured hexa
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
vol = G.getVolumeMap(a)
test.testT(vol, 4)

# 2D unstructured tri
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
vol = G.getVolumeMap(a)
test.testT(vol, 5)

# 3D unstructured tetra
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
vol = G.getVolumeMap(a)
test.testT(vol, 6)

# 3D unstructured pyra
a = G.cartPyra((0., 0., 0.), (1., 1., 1.), (5, 5, 5))
G._getVolumeMap(a)
test.testT(a, 9)

# 3D unstructured penta
a = G.cartPenta((0., 0., 0.), (1., 1., 1.), (5, 5, 5))
G._getVolumeMap(a)
test.testT(a, 10)

# 1D unstructured ME
a = G.cartHexa((0.,0.,0.), (0.1,0.0,0.0), (5,1,1))
b = G.cartHexa((0.4,0.,0.), (0.0,0.1,0.0), (1,11,1))
c = G.cartHexa((0.4,1.,0.), (0.1,0.0,0.), (4,1,1))
d = G.cartHexa((0.7,1.,0.), (0.1,0.1,0.05), (1,1,12))
a = C.mergeConnectivity([a, b, c, d], None)
G._getVolumeMap(a)
test.testT(a, 11)

# --- ME ---
# 2D unstructured ME: tri-quad-tri
#                          |
#                         tri
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (5,10,1))
b = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.2), (5,10,1))
c = G.cartTetra((0.8,0.,0.), (0.1,0.1,0.2), (5,10,1))
d = G.cartTetra((0.4,-0.9,0.), (0.1,0.1,0.2), (5,10,1))
a = C.mergeConnectivity([a, b, c, d], None)
G._getVolumeMap(a)
test.testT(a, 12)

# 3D unstructured ME: tetra   pyra
#                              |
#                     penta - hexa
a = G.cartTetra((-0.05,0.45,0.05), (0.1,0.1,0.1), (5,5,5))
b = G.cartPyra((0.4,0.4,0.), (0.1,0.1,0.1), (5,5,5))
c = G.cartPenta((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
d = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.mergeConnectivity([a, b, c, d], None)
G._getVolumeMap(a)
test.testT(a, 13)
