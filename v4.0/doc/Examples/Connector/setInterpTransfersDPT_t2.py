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
m = G.cartNGon((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
C._initVars(m, '{F}={CoordinateX}+2*{CoordinateY}+3*{CoordinateZ}')
C._initVars(m, '{G}=10*{CoordinateX}')
C._initVars(m,'centers:G',1.)
# Receiver mesh
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk)); a[0] = 'extraction'
C._initVars(a, 'F',-1000.)
C._initVars(a, 'G',-1000.)

# 2nd order, nodes, direct storage
t = C.newPyTree(['Rcv','Dnr']); t[2][1][2] = [a]; t[2][2][2] = [m]
C._initVars(t[2][1], 'cellN', 2)
t[2][2]=X.setInterpData(t[2][1],t[2][2],loc='nodes',storage='inverse',order=2,method='leastsquares')
info=X.setInterpTransfersD(t[2][2])
test.testO(info)
test.testA([info[0][1]],2)
