import Generator.PyTree as G
import Converter.PyTree as C
import Transform.PyTree as T
import Geom.PyTree as D

N = 4

# Donor mesh
a = G.cart((-2,-2,-2), (0.04,0.04,0.04), (200,100,100))
a = C.initVars(a, '{Density} = {CoordinateX}')
a = T.splitNParts(a, N=3*N)
t = C.newPyTree(['Base']) ; t[2][1][2] += a
C.convertPyTree2File(t, 'donor.cgns') 

# Receiver mesh
b = D.sphere((1,1,1), 0.5, N=100)
b = T.splitNParts(b, N=N)
t = C.newPyTree(['Base', b])
C.convertPyTree2File(t, 'receiver.cgns') 
