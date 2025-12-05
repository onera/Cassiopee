# - getVolumeMap (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# 3D structured
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,3))
a = C.addBC2Zone(a, 'wall1','BCWall','jmin')
a = C.addBC2Zone(a, 'match1','BCMatch','imin',a,'imax',[1,2,3])
a = C.addBC2Zone(a, 'match2','BCMatch','imax',a,'imin',[1,2,3])
a = C.fillEmptyBCWith(a,'overlap','BCOverlap')
t = C.newPyTree(['Base']); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.initVars(t,'Density',2.); t = C.initVars(t,'centers:cellN',1.)
t = G.getVolumeMap(t)
test.testT(t,1)

# 2D structured
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
t = C.newPyTree(['Base',2]); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.initVars(t,'Density',2.); t = C.initVars(t,'centers:cellN',1.)
t = G.getVolumeMap(t)
test.testT(t,2)

# 2D unstructured quad
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
t = C.newPyTree(['Base',2]); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.initVars(t,'Density',2.); t = C.initVars(t,'centers:cellN',1.)
t = G.getVolumeMap(t)
test.testT(t,3)

# 3D unstructured hexa
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
t = C.newPyTree(['Base']); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.initVars(t,'Density',2.); t = C.initVars(t,'centers:cellN',1.)
t = G.getVolumeMap(t)
test.testT(t,4)

# 2D unstructured tri
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
t = C.newPyTree(['Base',2]); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.initVars(t,'Density',2.); t = C.initVars(t,'centers:cellN',1.)
t = G.getVolumeMap(t)
test.testT(t,5)

# 3D unstructured tetra
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
t = C.newPyTree(['Base']); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.initVars(t,'Density',2.); t = C.initVars(t,'centers:cellN',1.)
t = G.getVolumeMap(t)
test.testT(t,6)
