# - XcellN (array) -
import Converter as C
import Connector as X
import Generator as G
import Geom as D
import Intersector as XOR
import KCore.test as test
import sys

# Test 1
# Mask
D = G.cart((0.,0.,0.), (0.5,0.4,0.6), (10,10,10))
D = C.convertArray2NGon(D)
#C.convertArrays2File([D], "a.plt")

# Mesh to blank
R = G.cart((-3.,-3.,-3.), (0.3,0.5,0.4), (20,20,20))
R = C.convertArray2NGon(R)
#C.convertArrays2File([R], "b.plt")

# celln init
cR = C.node2Center(R)
cR = C.initVars(cR, 'cellN', 2)
cR = C.extractVars(cR, ['cellN'])
#print cR

# Blanking
coef_and_indices = XOR.P1ConservativeChimeraCoeffs(R, cR, D)
test.testO(coef_and_indices,1)
