# - integNorm (array) -
import Converter as C
import Generator as G
import Post as P

ni = 30; nj = 40

# Maillage en noeuds
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))

# Champ a integrer en centres
c1 = C.array('vx,vy,vz', ni-1, nj-1, 1)
c = C.initVars(c1, 'vx,vy,vz', 1.); del c1

# Integration de chaque champ
res = P.integNorm([m], [c], [])
print(res)

# Champ a integrer en noeuds
c1 = C.array('vx,vy,vz', ni, nj, 1)
cn = C.initVars(c1, 'vx,vy,vz', 1.); del c1
resn = P.integNorm([m], [cn], [])
print(resn)
