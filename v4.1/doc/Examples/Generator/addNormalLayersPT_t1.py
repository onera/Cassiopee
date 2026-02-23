# - addNormalLayers (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

# sur une zone structuree
d = G.cart((0.1,0.,0.), (0.1,1,1),(2,1,1))
a = D.sphere((0,0,0), 1, 50 )
a = G.addNormalLayers(a, d)
test.testT(a,1)

# sur une surface TRI
a = D.sphere((0,0,0), 1, 50 )
a = C.convertArray2Tetra(a)
a = G.addNormalLayers(a, d)
test.testT(a,2)

# sur une surface QUAD
a = D.sphere((0,0,0), 1, 50 )
a = C.convertArray2Hexa(a)
a = G.addNormalLayers(a, d)
test.testT(a,3)

# sur une liste de zones structurees
d = G.cart((0.,0.,0.), (0.1,1,1),(3,1,1))
a = D.sphere6((0,0,0), 1, 20 )
a = C.initVars(a, 'Density',2.); a = C.initVars(a, 'centers:cellN',1.)
a = G.addNormalLayers(a, d)
test.testT(a,4)

# sur une liste de TRI
a = D.sphere6((0,0,0), 1, 20 ); a = C.convertArray2Tetra(a)
a = C.initVars(a, 'Density',2.); a = C.initVars(a, 'centers:cellN',1.)
a = G.addNormalLayers(a, d)
test.testT(a,5)

# sur une liste de QUAD
a = D.sphere6((0,0,0), 1, 20 ); a = C.convertArray2Hexa(a)
a = C.initVars(a, 'Density',2.); a = C.initVars(a, 'centers:cellN',1.)
a = G.addNormalLayers(a, d)
test.testT(a,6)

# sur un arbre
a = D.sphere6((0,0,0), 1, 20 )
t = C.newPyTree(['Base',2]); t[2][1][2] += a
t = G.addNormalLayers(t, d)
test.testT(a,7)
