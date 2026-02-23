# - extractMesh (array) -
import Converter as C
import Post as P
import Generator as G
import KCore.test as test

# Maillage en noeuds
ni = 30; nj = 40; nk = 10
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m, 'ro', 1.)

# Cree un maillage d'extraction
a = G.cart((0.,0.,0.), (0.01, 0.01, 0.1), (40, 40, 1))

# sur une liste - ordre 2 seulement
a2 = P.extractMesh([m], [a], order=2)
test.testA(a2,1)

# Extrait la solution sur le maillage d'extraction
for i in [2,3,5]:
    #print('Computing order %d...'%i)
    a2 = P.extractMesh([m], a, i)
    test.testA([m,a2], i)

# With a hook
hook = C.createHook([m], function='extractMesh')
for i in [2,3,5]:
    a2 = P.extractMesh([m], a, i, hook=[hook])
    test.testA([m,a2], i)
C.freeHook(hook)
