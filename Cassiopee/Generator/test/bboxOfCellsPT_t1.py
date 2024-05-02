# - bboxOfCells (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# test 1D structure
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(20,1,1))
a = C.initVars(a,'Density',2.); a = C.initVars(a,'centers:cellN',1.)
a = G.bboxOfCells(a)
test.testT(a,1)

# test 2D structure
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(20,20,1))
a = C.addBC2Zone(a, 'wall1','BCWall','jmin')
a = C.initVars(a,'Density',2.); a = C.initVars(a,'centers:cellN',1.)
a = G.bboxOfCells(a)
test.testT(a,2)

# test 3d structure
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(20,20,20))
a = C.addBC2Zone(a, 'wall1','BCWall','jmin')
a = C.addBC2Zone(a, 'match1','BCMatch','imin',a,'imax',[1,2,3])
a = C.addBC2Zone(a, 'match2','BCMatch','imax',a,'imin',[1,2,3])
a = C.fillEmptyBCWith(a,'overlap','BCOverlap')
a = C.initVars(a,'Density',2.); a = C.initVars(a,'centers:cellN',1.)
a = G.bboxOfCells(a)
test.testT(a,3)

# test TRI
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(20,20,1))
a = C.convertArray2Tetra(a); a = G.bboxOfCells(a)
a = C.initVars(a,'Density',2.); a = C.initVars(a,'centers:cellN',1.)
test.testT(a,4)

# test QUAD
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(20,20,1))
a = C.convertArray2Hexa(a); a = G.bboxOfCells(a)
a = C.initVars(a,'Density',2.); a = C.initVars(a,'centers:cellN',1.)
test.testT(a,5)

# test TETRA
a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))
a = C.convertArray2Tetra(a); a = G.bboxOfCells(a)
a = C.initVars(a,'Density',2.); a = C.initVars(a,'centers:cellN',1.)
test.testT(a,6)

# test HEXA
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(10,10,10))
a = C.convertArray2Hexa(a); a = G.bboxOfCells(a)
a = C.initVars(a,'Density',2.); a = C.initVars(a,'centers:cellN',1.)
test.testT(a,7)
