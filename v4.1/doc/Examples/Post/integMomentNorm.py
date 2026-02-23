# - integMomentNorm (array) -
import Converter as C
import Generator as G
import Post as P

# Maillage et champs non structure, en noeuds
m = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
c = C.array('ro', 100, 162, 'TRI')
c = C.initVars(c, 'ro', 1.)
res = P.integMomentNorm([m], [c], [], (5.,5., 0.)); print(res)

# Maillage en noeuds
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))

# Champ a integrer en centres
c = C.array('v', ni-1, nj-1, 1)
c = C.initVars(c, 'v', 1.)

# Integration de chaque champ
res = P.integMomentNorm([m], [c], [], (5.,5., 0.) ); print(res)

# Champ a integrer en noeuds
cn = C.array('v', ni, nj, 1)
cn = C.initVars(cn, 'v', 1.)
resn = P.integMomentNorm([m], [cn], [], (5.,5., 0.)); print(resn)
