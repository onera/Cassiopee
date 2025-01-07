# - integ (array) -
import Converter as C
import Generator as G
import Post as P

# Node mesh
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))

# Field in centers
c = C.array('vx', ni-1, nj-1, 1); c = C.initVars(c, 'vx', 1.)
resc = P.integ([m], [c], [])[0]; print(resc)

# Field in nodes
cn = C.array('vx', ni, nj, 1); cn = C.initVars(cn, 'vx', 1.)
resn = P.integ([m], [cn], [])[0]; print(resn)
