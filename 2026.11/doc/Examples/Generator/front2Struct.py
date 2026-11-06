# - front2Struct (array) -
import Geom as D
import Generator as G
import Converter as C
import Connector as X
import Post as P
import Transform as T

# Surface
a = D.sphere( (0,0,0), 1, 5 ); a = C.convertArray2Tetra(a)

# Grille cartesienne
ni = 10; nj = 10; nk = 10
hi = 2.2/(ni-1); hj = 2.2/(nj-1); hk = 2.2/(nk-1)
b = G.cart( (-1.-0.1, -1.-0.1, -1.-0.1), (hi, hj, hk), (ni,nj,nk) )
celln = C.array('cellN', ni, nj, nk)
celln = C.initVars(celln, 'cellN', 1.)

# Masquage
cellno = X.blankCells([b], [celln], [a], blankingType=0, delta=1.e-10, dim=3)
a = C.initVars(a, 'cellN', 1)
b = C.addVars([b, cellno[0]])

# Selection du front
b1 = P.selectCells2(b, cellno[0])
bext = P.exteriorFaces(b1)
bexts = T.splitConnexity(bext)
f = bexts[1]; f = C.initVars(f, 'cellN', 1)

# Lissage du front
f = T.smooth(f, eps=0.5, niter=5) # lissage du front quad
#C.convertArrays2File([f], 'out.plt')

# Creation du maillage
distrib = G.cart( (0,0,0), (0.01,1,1), (10,1,1) )
m = G.front2Struct(f, a, distrib, 10)

C.convertArrays2File(m, 'out.plt')
