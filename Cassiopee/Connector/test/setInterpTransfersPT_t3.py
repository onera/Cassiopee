# - setInterpTransfers (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import KCore.test as test
# DONOR NON STRUCTURE TETRA
# Create a function
def F(x,y,z):
    deg = 1
    if deg == 0: return 10.
    elif deg == 1: return x + 2.*y + 3.*z
    elif deg == 2: return x*x + 2.*y*y + 3*z
    elif deg == 3: return x*x*y + 2.*y*y*y + 3*z
    elif deg == 4: return x*x*x*x + 2.*y*y*y*y +z*z
    else: return 2*x*x*x*x*x + 2.*y*y*z + z*z

# Donor mesh
ni = 11; nj = 11; nk = 11
m = G.cartTetra((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
C._initVars(m, 'F', F, ['CoordinateX','CoordinateY','CoordinateZ'])
C._initVars(m,'centers:G',1.)

# Receiver mesh
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk)); a[0] = 'extraction'

# 2nd order, nodes, direct storage
t = C.newPyTree(['Rcv',a,'Dnr',m])
t[2][1]=C.initVars(t[2][1], 'cellN', 2) # tag pour interpoler tous les pts
t[2][1]=X.setInterpData(t[2][1],t[2][2],loc='nodes',storage='direct',order=2)
t[2][1]=C.initVars(t[2][1],'F',0.)
t[2][1]=X.setInterpTransfers(t[2][1],t[2][2],variables=['F'])
C._rmVars(t[2][1], 'cellN')
test.testT(t,1)
#
t = C.newPyTree(['Rcv',a,'Dnr',m])
C._initVars(t[2][1], 'centers:cellN', 2) # tag pour interpoler tous les pts
t[2][1]=X.setInterpData(t[2][1],t[2][2],loc='centers',storage='direct',order=2)
C._rmVars(t[2][1],['F'])
C._initVars(t[2][1],'centers:F',0.)
X._setInterpTransfers(t[2][1],t[2][2],variables=['F'])
C._rmVars(t[2][1], 'centers:cellN')
test.testT(t,2)
