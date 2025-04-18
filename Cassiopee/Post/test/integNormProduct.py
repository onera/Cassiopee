# - integNormProduct (array) -
import Converter as C
import Generator as G
import Post as P

# Maillage et champs non structure, en noeuds
m = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
c = C.array('vx,vy,vz', m[1].shape[1], m[2].shape[1], 'TRI')
c = C.initVars(c, 'vx,vy,vz', 1.)
res = P.integNormProduct([m], [c], []); print(res)

# Maillage en noeuds
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))

# Champ a integrer en centres
c = C.array('vx,vy,vz', ni-1, nj-1, 1)
c = C.initVars(c, 'vx,vy,vz', 1.)

# Integration de chaque champ
res = P.integNormProduct([m], [c], []); print(res)

# Champ a integrer en noeuds
cn = C.array('vx,vy,vz', ni, nj, 1)
cn = C.initVars(cn, 'vx,vy,vz', 1.)
resn = P.integNormProduct([m], [cn], []); print(resn)
