# - integ (array) -
import Converter as C
import Generator as G
import Post as P
import KCore.test as test

ni = 11; nj = 11

def f1(x,y):
    return 2*x + y

def f2(x,y) :
    return 3*x*y + 4

# STRUCT 2D NODE / CENTER
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m2 = C.node2Center(m)
c = C.array('vx,vy,vz', ni-1, nj-1, 1)
c0 = C.addVars([m2,c])
c0 = C.initVars(c0,'vx', f1, ['x','y'])
c0 = C.initVars(c0,'vy', f2, ['x','y'])
c = C.extractVars(c0, ['vx','vy','vz'])
res = P.integ([m], [c], [])
out = C.array('res', 3, 1, 1)
out[1][0,0] = res[0]; out[1][0,1] = res[1]; out[1][0,2] = res[2];
test.testA([out], 1)

# STRUCT 2D NODE / NODE
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
c = C.array('vx,vy,vz', ni, nj, 1)
c0 = C.addVars([m,c])
c0 = C.initVars(c0,'vx', f1, ['x','y'])
c0 = C.initVars(c0,'vy', f2, ['x','y'])
c = C.extractVars(c0, ['vx','vy','vz'])
res = P.integ([m], [c], [])
out = C.array('res', 3, 1, 1)
out[1][0,0] = res[0]; out[1][0,1] = res[1]; out[1][0,2] = res[2];
test.testA([out], 2)

# STRUCT 1D NODE / CENTER
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,1,1))
m2 = G.cart((0.5,0,0), (9./(ni-2),9./(nj-2),1), (ni-1,1,1))
c = C.array('vx,vy,vz', ni-1, 1, 1)
c0 = C.addVars([m2,c])
c0 = C.initVars(c0,'vx', f1, ['x','y'])
c0 = C.initVars(c0,'vy', f2, ['x','y'])
c = C.extractVars(c0, ['vx','vy','vz'])
res = P.integ([m], [c], [])
out = C.array('res', 3, 1, 1)
out[1][0,0] = res[0]; out[1][0,1] = res[1]; out[1][0,2] = res[2];
test.testA([out], 3)

# STRUCT 1D NODE / NODE
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,1,1))
c = C.array('vx,vy,vz', ni, 1, 1)
c0 = C.addVars([m,c])
c0 = C.initVars(c0,'vx', f1, ['x','y'])
c0 = C.initVars(c0,'vy', f2, ['x','y'])
c = C.extractVars(c0, ['vx','vy','vz'])
res = P.integ([m], [c], [])
out = C.array('res', 3, 1, 1)
out[1][0,0] = res[0]; out[1][0,1] = res[1]; out[1][0,2] = res[2];
test.testA([out], 4)

# TRI NODE / NODE
m = G.cartTetra((0.,0.,0.),(10./(ni-1),10./(nj-1),1), (ni,nj,1))
c = C.array('vx,vy,vz', m[1].shape[1], m[2].shape[1], 'TRI')
c = C.initVars(c, 'vx,vy,vz', 1.)
res = P.integ([m], [c], [])
out = C.array('res', 3, 1, 1)
out[1][0,0] = res[0]; out[1][0,1] = res[1]; out[1][0,2] = res[2];
test.testA([out], 5)

# TRI NODE / CENTER
m = G.cartTetra((0.,0.,0.), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
c = C.array('vx,vy,vz', m[1].shape[1], m[2].shape[1], 'TRI')
c = C.initVars(c, 'vx,vy,vz', 1.)
c = C.node2Center(c)
res = P.integ([m], [c], [])
out = C.array('res', 3, 1, 1)
out[1][0,0] = res[0]; out[1][0,1] = res[1]; out[1][0,2] = res[2];
test.testA([out], 6)

# BAR NODE / NODE
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,1,1))
m = C.convertArray2Tetra(a)
c = C.array('vx,vy,vz', m[1].shape[1], m[2].shape[1], 'BAR')
c = C.initVars(c, 'vx,vy,vz', 1.)
res = P.integ([m], [c], [])
out = C.array('res', 3, 1, 1)
out[1][0,0] = res[0]; out[1][0,1] = res[1]; out[1][0,2] = res[2];
test.testA([out], 7)

# BAR NODE / CENTER
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,1,1))
m = C.convertArray2Tetra(a)
c = C.array('vx,vy,vz', m[1].shape[1], m[2].shape[1], 'BAR')
c = C.initVars(c, 'vx,vy,vz', 1.)
c = C.addVars([m,c])
c = C.node2Center(c)
c = C.extractVars(c,['vx','vy','vz'])
res = P.integ([m], [c], [])
out = C.array('res', 3, 1, 1)
out[1][0,0] = res[0]; out[1][0,1] = res[1]; out[1][0,2] = res[2];
test.testA([out], 8)
