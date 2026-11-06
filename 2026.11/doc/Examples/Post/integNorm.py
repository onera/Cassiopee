# - integNorm (array) -
import Converter as C
import Generator as G
import Post as P

# Node mesh and field
m = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
c1 = C.array('ro', 100, 162, 'TRI')
c = C.initVars(c1, 'ro', 1.)
res = P.integNorm([m], [c], []); print('res1 = %f'%res)

# Node mesh
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))

# Centers field
c1 = C.array('vx', ni-1, nj-1, 1)
c = C.initVars(c1, 'vx', 1.)

# Integration
res = P.integNorm([m], [c], []); print('res2 = %f'%res)

# Node field
c1 = C.array('vx, vy', ni, nj, 1)
cn = C.initVars(c1, 'vx', 1.)
cn = C.initVars(c1, 'vy', 1.)
resn = P.integNorm([m], [cn], []); print('res3 = %f'%resn)
