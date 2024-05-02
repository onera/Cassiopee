# - integMomentNorm -
import Converter as C
import Generator as G
import Post as P
import KCore.test as test

ni = 11; nj = 11
def f1(x,y):
    return 2*x + y

def f2(x,y) :
    return 3*x*y + 4

xc = 5.; yc = 1.; zc = 0.

# STRUCT 2D NODE / CENTER
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m2 = C.node2Center(m)
c = C.array('vx,vy', ni-1, nj-1, 1)
c0 = C.addVars([m2,c])
c0 = C.initVars(c0,'vx', f1, ['x','y'])
c0 = C.initVars(c0,'vy', f2, ['x','y'])
c = C.extractVars(c0, ['vx','vy'])
res = P.integMomentNorm([m], [c], [], (xc,yc,zc))
out = C.array('resvx,resvy', 3, 1, 1)
out[1][0,0] = res[0][0]; out[1][0,1] = res[0][1]; out[1][0,2] = res[0][2];
out[1][0,0] = res[1][0]; out[1][0,1] = res[1][1]; out[1][0,2] = res[1][2];
test.testA([out], 1)

# STRUCT 2D NODE / NODE
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
c = C.array('vx,vy', ni-1, nj-1, 1)
c0 = C.addVars([m2,c])
c0 = C.initVars(c0,'vx', f1, ['x','y'])
c0 = C.initVars(c0,'vy', f2, ['x','y'])
c = C.extractVars(c0, ['vx','vy'])
res = P.integMomentNorm([m], [c], [], (xc,yc,zc))
out = C.array('resvx,resvy', 3, 1, 1)
out[1][0,0] = res[0][0]; out[1][0,1] = res[0][1]; out[1][0,2] = res[0][2];
out[1][0,0] = res[1][0]; out[1][0,1] = res[1][1]; out[1][0,2] = res[1][2];
test.testA([out], 2)


# TRI NODE / NODE
m = G.cartTetra((0.,0.,0.),(10./(ni-1),10./(nj-1),1), (ni,nj,1))
c = C.array('vx,vy', m[1].shape[1], m[2].shape[1], 'TRI')
c0 = C.addVars([m,c])
c0 = C.initVars(c0,'vx', 1.)
c0 = C.initVars(c0,'vy', 1.)
c = C.extractVars(c0, ['vx','vy'])
res = P.integMomentNorm([m], [c], [], (xc,yc,zc))
out = C.array('resvx,resvy', 3, 1, 1)
out[1][0,0] = res[0][0]; out[1][0,1] = res[0][1]; out[1][0,2] = res[0][2];
out[1][0,0] = res[1][0]; out[1][0,1] = res[1][1]; out[1][0,2] = res[1][2];
test.testA([out], 3)

# TRI NODE / CENTER
c = C.array('vx,vy', m[1].shape[1], m[2].shape[1], 'TRI')
c0 = C.addVars([m,c])
c0 = C.initVars(c0,'vx', 1.)
c0 = C.initVars(c0,'vy', 1.)
c0 = C.node2Center(c0)
c = C.extractVars(c0, ['vx','vy'])

res = P.integMomentNorm([m], [c], [], (xc,yc,zc))
out = C.array('resvx,resvy', 3, 1, 1)
out[1][0,0] = res[0][0]; out[1][0,1] = res[0][1]; out[1][0,2] = res[0][2];
out[1][0,0] = res[1][0]; out[1][0,1] = res[1][1]; out[1][0,2] = res[1][2];
test.testA([out], 4)
