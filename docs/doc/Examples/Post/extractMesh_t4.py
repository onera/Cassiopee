# - extractMesh
#  maillage d extraction non structure -
import Converter as C
import Post as P
import Generator as G
import KCore.test as test

# Maillage en noeuds
ni = 30; nj = 40; nk = 10
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))

c = C.array('ro,rou, rov,row,roe,cellN', ni, nj, nk)
c = C.initVars(c, 'ro', 1.)
c = C.initVars(c, 'rou', 1.)
c = C.initVars(c, 'rov', 0.)
c = C.initVars(c, 'row', 0.)
c = C.initVars(c, 'roe', 1.)
c = C.initVars(c, 'cellN', 1)
m = C.addVars([m,c])

# Cree un maillage d'extraction non structure
a = G.cart( (0.,0.,0.), (1., 0.1, 0.1), (20, 20, 1))
a = C.convertArray2Tetra(a)

# sur une liste - ordre 2 seulement
a2 = P.extractMesh([m], [a])
test.testA(a2,1)

# Extrait la solution sur le maillage d'extraction
for i in [2,3,5]:
    print('Computing order %d...'%i)
    a2 = P.extractMesh([m], a, i)
    test.testA([a2], i)
