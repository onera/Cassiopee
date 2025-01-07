# - close (pyTree)-
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test
import Geom.PyTree as D
import Transform.PyTree as T

# test 1D : circle
a = D.circle((0.,0.,0.), 1., 0., 359.995, 500)
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
a = C.convertArray2Tetra(a)
a2 = G.close(a, 1.e-4)
test.testT(a2,1)

# test 2D cylindre QUAD
ni = 20; nj = 20; nk = 5
a0 = G.cylinder((0.,0.,0.), 0., 1., 0., 359, 1., (ni,nj,nk))
a = T.subzone(a0, (1,nj,1),(ni,nj,nk))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
a = C.convertArray2Hexa(a)
a2 = G.close(a, 1.e-1)
test.testT(a2,2)

# test 2D TRI
a = T.subzone(a0, (1,nj,1),(ni,nj,nk))
a = C.convertArray2Tetra(a)
C._addVars(a,'Density'); C._initVars(a,'centers:cellN',1.)
a2 = G.close(a, 1.e-1)
test.testT(a2,3)

# test 3D cylindre HEXA
ni = 20; nj = 20; nk = 5
a0 = G.cylinder((0.,0.,0.), 0., 1., 0., 359, 1., (ni,nj,nk))
a0 = T.subzone(a0,(1,10,1),(20,13,5))
a = C.convertArray2Hexa(a0)
C._addVars(a,'Density'); C._initVars(a,'centers:cellN',1.)
a2 = G.close(a, 1.e-2)
test.testT(a2,4)

# test 3D TETRA
a = C.convertArray2Tetra(a0)
C._addVars(a,'Density'); C._initVars(a,'centers:cellN',1.)
a2 = G.close(a,2.e-2)
test.testT(a2,5)

# test 3D PENTA
a = G.cartPenta((0.,0.,0.), (0.1,0.1,0.2), (10,10,10)); a = C.addVars(a,'F')
C._addVars(a,'Density'); C._initVars(a,'centers:cellN',1.)
a2 = G.close(a, 2.e-2)
test.testT([a2],6)

# test 3D PYRA
a = G.cartPyra((0.,0.,0.), (0.1,0.1,0.2), (10,10,10)); a = C.addVars(a,'F')
C._addVars(a,'Density'); C._initVars(a,'centers:cellN',1.)
a2 = G.close(a, 2.e-2)
test.testT([a2],7)
