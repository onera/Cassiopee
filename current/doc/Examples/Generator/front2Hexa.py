# - front2Hexa (array) -
import Geom as D
import Generator as G
import Converter as C
import Connector as X
import Post as P
import Transform as T

# Surface
a = D.sphere((0,0,0), 1); a = C.convertArray2Tetra(a)

# Grille cartesienne
ni = 30; nj = 30; nk = 30
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
z = bexts[1]
z = C.initVars(z, 'cellN', 1)

# Projection du front
t = G.front2Hexa(z, a, 0.01, 0.01, 0.1, 50)

C.convertArrays2File([t], 'out.plt')
