# - snapFront (array) -
import Generator as G
import Converter as C
import Geom as D
import Connector as X
import Transform as T

s = D.circle((0,0,0), 1., N=100)
s = T.addkplane(s)

# Grille cartesienne (reguliere)
BB = G.bbox([s])
ni = 100; nj = 100; nk = 3
xmin = BB[0]; ymin = BB[1]; zmin = BB[2]-0.5
xmax = BB[3]; ymax = BB[4]; zmax = BB[5]+0.5
hi = (xmax-xmin)/(ni-1); hj = (ymax-ymin)/(nj-1)
h = min(hi, hj)
ni = int((xmax-xmin)/h)+7; nj = int((ymax-ymin)/h)+7
b = G.cart((xmin-3*h, ymin-3*h, zmin), (h, h, 1.), (ni,nj,nk))
celln = C.array('cellN', ni, nj, nk)
celln = C.initVars(celln, 'cellN', 1.)

# Masquage
cellno = X.blankCells([b], [celln], [s], blankingType=0, delta=0., dim=2)
a = C.initVars(s, 'cellN', 1)
b = C.addVars([b, cellno[0]])

# Adapte le front de la grille a la surface
b = T.subzone(b, (1,1,2), (b[2],b[3],2))
b = G.snapFront(b, [s])

C.convertArrays2File([a,b], 'out.plt')
