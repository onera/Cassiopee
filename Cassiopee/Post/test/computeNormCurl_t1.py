# - computeNormCurl (array) -
import Converter as C
import Post as P
import Generator as G
import KCore.test as T

def F(x,y,z):
    return 12*y*y + 4

def DF(y):
    return -24.*y

# Test 3D
ni = 30; nj = 40; nk = 3
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m,'u',F,['x','y','z'])
m = C.initVars(m,'v',0.)
m = C.initVars(m,'w',0.)
m1 = C.initVars(m,'sol',DF,['y'])
p1 = P.computeNormCurl(m1,['u','v','w']) # defined on centers
T.testA([p1], 1)

# Test 2D
ni = 30; nj = 40; nk = 1
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m,'u',F,['x','y','z'])
m = C.initVars(m,'v',0.)
m2 = C.initVars(m,'w',0.)
p2 = P.computeNormCurl(m2,['u','v','w']) # defined on centers
T.testA([p2], 2)

# test sur une liste
P = P.computeNormCurl([m1,m2],['u','v','w']) # defined on centers
T.testA(P, 3)
