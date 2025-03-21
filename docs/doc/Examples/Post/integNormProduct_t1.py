# - integNormProduct (array) -
import Converter as C
import Generator as G
import Post as P
import KCore.test as test

ni = 11; nj = 11

def f3(z) :
    return 2*z

# STRUCT 2D NODE / CENTER
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m2 = C.node2Center(m)
c = C.array('vx,vy,vz', ni-1, nj-1, 1)
c0 = C.addVars([m2,c])
c0 = C.initVars(c0,'vz', 1.)
c = C.extractVars(c0, ['vx','vy','vz'])
res = P.integNormProduct([m], [c], [])
out = C.array('res', 1, 1, 1)
out[1][0][0] = res
test.testA([out], 1)

# STRUCT 2D NODE / NODE
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
c = C.array('vx,vy,vz', ni, nj, 1)
c0 = C.addVars([m,c])
c0 = C.initVars(c0,'vz', 1.)
c = C.extractVars(c0, ['vx','vy','vz'])
res = P.integNormProduct([m], [c], [])
out = C.array('res', 1, 1, 1)
out[1][0][0] = res
test.testA([out], 2)

# TRI NODE / NODE
m = G.cartTetra((0.,0.,0.),(10./(ni-1),10./(nj-1),1), (ni,nj,1))
c = C.array('vx,vy,vz', m[1].shape[1], m[2].shape[1], 'TRI')
c0 = C.addVars([m,c])
c0 = C.initVars(c0,'vz', 1.)
c = C.extractVars(c0, ['vx','vy','vz'])
res = P.integNormProduct([m], [c], [])
out = C.array('res', 1, 1, 1)
out[1][0][0] = res
test.testA([out], 3)

# TRI NODE / CENTER
m = G.cartTetra((0.,0.,0.),(10./(ni-1),10./(nj-1),1), (ni,nj,1))
c = C.array('vx,vy,vz', m[1].shape[1], m[2].shape[1], 'TRI')
c0 = C.addVars([m,c])
c0 = C.initVars(c0,'vz', 1.)
c0 = C.node2Center(c0)
c = C.extractVars(c0, ['vx','vy','vz'])
res = P.integNormProduct([m], [c], [])
out = C.array('res', 1, 1, 1)
out[1][0][0] = res
test.testA([out], 4)
