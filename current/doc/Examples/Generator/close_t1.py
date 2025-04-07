# - close (array) -
import Converter as C
import Generator as G
import KCore.test as test

a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0.01, 10., (20,20,10))
a = C.addVars(a, 'F')

# test structure
a1 = G.close(a, 1.e-3)
test.testA([a1], 1)

# test non structure hexa
a2 = C.convertArray2Hexa(a)
a2 = G.close(a2, 1.e-3)
test.testA([a2], 2)

# test close non structure tetra
a3 = C.convertArray2Tetra(a)
a3 = G.close(a3, 1.e-3)
test.testA([a3], 3)

# test close NGON avec retour de la table d indir. des vertices
indices = []
a4 =  G.cylinder((0.,0.,0.), 0.5, 1., 360., 0.01, 10., (10,10,10))
a4 = C.convertArray2NGon(a4)
a4 = C.addVars(a4, 'F')
a4 = G.close(a4, 5.e-3, indices=indices)
test.testA([a4], 4)
test.testO(indices, 5)
