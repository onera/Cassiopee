# - selectInsideElts (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,1))
a = C.convertArray2Tetra(a)
C._addVars(a,'Density'); C._initVars(a,'centers:cellN',1)
b = D.circle((5,5,0), 3.)
b = C.convertArray2Tetra(b)
a2 = G.selectInsideElts(a, b); a2[0] = 'cart2'
test.testT(a)
