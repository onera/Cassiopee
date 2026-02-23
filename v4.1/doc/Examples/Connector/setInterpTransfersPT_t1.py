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

# 2nd order, nodes, direct storage
t = C.newPyTree(['Rcv','Dnr']); t[2][1][2] = [a]; t[2][2][2] = [m]
t[2][1] = C.initVars(t[2][1], 'cellN', 2)
t[2][1] = X.setInterpData(t[2][1],t[2][2],loc='nodes',storage='direct',order=2)
t[2][1] = C.initVars(t[2][1],'F',0.)
t[2][1] = X.setInterpTransfers(t[2][1],t[2][2],variables=['F'])
t[2][1] = C.rmVars(t[2][1], 'cellN')
test.testT(t,1)

# 3rd order, nodes, direct storage
t = C.newPyTree(['Rcv','Dnr']); t[2][1][2] = [a]; t[2][2][2] = [m]
t[2][1] = C.initVars(t[2][1], 'cellN', 2)
t[2][1] = X.setInterpData(t[2][1],t[2][2],loc='nodes',storage='direct',order=3)
t[2][1] = C.initVars(t[2][1],'F',0.)
t[2][1] = X.setInterpTransfers(t[2][1],t[2][2],variables=['F'])
t[2][1] = C.rmVars(t[2][1], 'cellN')
test.testT(t,2)

# 5th order direct storage, at nodes
t = C.newPyTree(['Rcv','Dnr']); t[2][1][2] = [a]; t[2][2][2] = [m]
t[2][1] = C.initVars(t[2][1], 'cellN', 2)
t[2][1] = X.setInterpData(t[2][1],t[2][2],loc='nodes',storage='direct',order=5)
t[2][1] = C.initVars(t[2][1],'F',0.)
t[2][1] = X.setInterpTransfers(t[2][1],t[2][2],variables=['F'])
t[2][1] = C.rmVars(t[2][1], 'cellN')
test.testT(t,3)
