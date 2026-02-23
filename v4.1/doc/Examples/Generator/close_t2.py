# - close (array) -
import Generator as G
import Converter as C
import KCore.test as test
import Geom as D
import Transform as T

# test 1D : circle
a = D.circle((0.,0.,0.), 1., 0., 359.995, 500)
a = C.convertArray2Tetra(a)
a = C.addVars(a,'F')
a2 = G.close(a, 1.e-4)
test.testA([a2],1)

# test 2D cylindre QUAD
a0 = G.cylinder((0.,0.,0.), 0., 1., 0., 359, 1., (20,20,5))
a0 = T.subzone(a0, (1,a0[3],1),(a0[2],a0[3],a0[4]))
a = C.convertArray2Hexa(a0); a = C.addVars(a,'F')
a2 = G.close(a, 0.1)
test.testA([a2],2)

# test 2D TRI
a = C.convertArray2Tetra(a0)
a = C.addVars(a,'F')
a2 = G.close(a, 0.1)
test.testA([a2],3)

# test 3D cylindre HEXA
a0 = G.cylinder((0.,0.,0.), 0., 1., 0., 359, 1., (20,20,5))
a0 = T.subzone(a0,(1,10,1),(20,13,5))
a = C.convertArray2Hexa(a0); a = C.addVars(a,'F')
a2 = G.close(a, 0.01)
test.testA([a2],4)

# test 3D TETRA
a = C.convertArray2Tetra(a0); a = C.addVars(a,'F')
a2 = G.close(a, 2.e-2)
test.testA([a2],5)

# test 3D PENTA
a = G.cartPenta((0.,0.,0.), (0.1,0.1,0.2), (10,10,10)); a = C.addVars(a,'F')
a2 = G.close(a, 2.e-2)
test.testA([a2],6)

# test 3D PYRA
a = G.cartPyra((0.,0.,0.), (0.1,0.1,0.2), (10,10,10)); a = C.addVars(a,'F')
a2 = G.close(a, 2.e-2)
test.testA([a2],7)
