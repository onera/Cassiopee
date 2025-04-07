# - setInterpTransfers (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import KCore.test as test

# Create a function
def F(x,y,z):
    deg = 1
    if deg == 0 : return 10.
    elif deg == 1 : return x + 2.*y + 3.*z
    elif deg == 2 : return x*x + 2.*y*y + 3*z
    elif deg == 3 : return x*x*y + 2.*y*y*y + 3*z
    elif deg == 4 : return x*x*x*x + 2.*y*y*y*y +z*z
    else : return 2*x*x*x*x*x + 2.*y*y*z + z*z

# Donor mesh
ni = 11; nj = 11; nk = 11
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m, 'F', F, ['CoordinateX','CoordinateY','CoordinateZ'])
m = C.initVars(m,'centers:G',1.)
# Receiver mesh
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk)); a[0] = 'extraction'

# 2nd order, centers, direct storage
t = C.newPyTree(['Rcv','Dnr']); t[2][1][2] = [a]; t[2][2][2] = [m]
t[2][1] = C.initVars(t[2][1], '{centers:cellN}=2')
t[2][1] = X.setInterpData(t[2][1],t[2][2],loc='centers',storage='direct')
t[2][1] = C.initVars(t[2][1],'centers:F',0.)
t[2][1] = X.setInterpTransfers(t[2][1],t[2][2],variables=['F'])
test.testT(t,1)

# 2nd order, centers, inverse storage
t = C.newPyTree(['Rcv','Dnr']); t[2][1][2] = [a]; t[2][2][2] = [m]
t[2][1] = C.initVars(t[2][1], '{centers:cellN}=2')
t[2][2] = X.setInterpData(t[2][1],t[2][2],loc='centers',storage='inverse')
t[2][1] = C.initVars(t[2][1],'centers:F',0.)
t[2][1] = X.setInterpTransfers(t[2][1],t[2][2],variables=['F'])
test.testT(t,2)

# 2nd order, nodes, direct storage, sans variables en entree
t = C.newPyTree(['Rcv','Dnr']); t[2][1][2] = [a]; t[2][2][2] = [m]
t[2][1] = C.initVars(t[2][1], 'cellN',2)
t[2][1] = X.setInterpData(t[2][1],t[2][2],loc='nodes',storage='direct')
t[2][1] = C.initVars(t[2][1],'F',0.)
t[2][1] = X.setInterpTransfers(t[2][1],t[2][2])
test.testT(t,3)
