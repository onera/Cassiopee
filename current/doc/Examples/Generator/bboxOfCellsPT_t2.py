# - bboxOfCells (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# test 1D structure
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(20,1,1))
t = C.newPyTree(['Base',1]); t[2][1][2].append(a)
t = C.initVars(t,'Density',2.); t = C.initVars(t,'centers:cellN',1.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = G.bboxOfCells(t)
test.testT(t,1)

# test 2D structure
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(20,20,1))
a = C.addBC2Zone(a, 'wall1','BCWall','jmin')
t = C.newPyTree(['Base',2]); t[2][1][2].append(a)
t = C.initVars(t,'Density',2.); t = C.initVars(t,'centers:cellN',1.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = G.bboxOfCells(t)
test.testT(t,2)

# test 3d structure
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(20,20,20))
a = C.addBC2Zone(a, 'wall1','BCWall','jmin')
a = C.addBC2Zone(a, 'match1','BCMatch','imin',a,'imax',[1,2,3])
a = C.addBC2Zone(a, 'match2','BCMatch','imax',a,'imin',[1,2,3])
a = C.fillEmptyBCWith(a,'overlap','BCOverlap')
t = C.newPyTree(['Base']); t[2][1][2].append(a)
t = C.initVars(t,'Density',2.); t = C.initVars(t,'centers:cellN',1.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = G.bboxOfCells(t)
test.testT(t,3)

# test TRI
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(20,20,1))
a = C.convertArray2Tetra(a);
t = C.newPyTree(['Base',2]); t[2][1][2].append(a)
t = C.initVars(t,'Density',2.); t = C.initVars(t,'centers:cellN',1.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = G.bboxOfCells(t)
test.testT(t,4)

# test QUAD
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(20,20,1))
a = C.convertArray2Hexa(a)
t = C.newPyTree(['Base',2]); t[2][1][2].append(a)
t = C.initVars(t,'Density',2.); t = C.initVars(t,'centers:cellN',1.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = G.bboxOfCells(t)
test.testT(t,5)

# test TETRA
a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))
a = C.convertArray2Tetra(a)
t = C.newPyTree(['Base']); t[2][1][2].append(a)
t = C.initVars(t,'Density',2.); t = C.initVars(t,'centers:cellN',1.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = G.bboxOfCells(t)
test.testT(t,6)

# test HEXA
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(10,10,10))
a = C.convertArray2Hexa(a)
t = C.newPyTree(['Base']); t[2][1][2].append(a)
t = C.initVars(t,'Density',2.); t = C.initVars(t,'centers:cellN',1.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = G.bboxOfCells(t)
test.testT(t,7)
