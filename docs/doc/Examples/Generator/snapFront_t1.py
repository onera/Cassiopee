# - snapFront (array) -
import Generator as G
import Converter as C
import Geom as D
import Connector as X
import Transform as T
import KCore.test as test
import Post as P

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
b = G.cart( (xmin-3*h, ymin-3*h, zmin), (h, h, 1.), (ni,nj,nk) )
celln = C.array('cellN', ni, nj, nk)
celln = C.initVars(celln, 'cellN', 1.)

# Masquage
cellno = X.blankCells([b], [celln], [s], blankingType=0, delta=0., dim=2)
a = C.initVars(s, 'cellN', 1)
b = C.addVars([b, cellno[0]])

# Adapte le front de la grille a la surface
b = T.subzone(b, (1,1,2), (b[2],b[3],2))
f = G.snapFront(b, [s])

test.testA([f], 1)

# Adapte le front de la grille a la surface avec optimisation du front
f = G.snapFront(b, [s], optimized=1)

test.testA([f], 11)

# Grille non structuree
b = C.convertArray2Hexa(b)
f = G.snapFront(b, [s])
test.testA([f], 2)
# Adapte le front de la grille a la surface avec optimisation du front
f = G.snapFront(b, [s], optimized=1)
test.testA([f], 21)

s = D.polyline([(0.02,0,0),(1,1,0),(2,1,0),(0.02,0,0)])
s1 = T.addkplane(s)
contours = P.exteriorFaces(s1)
contours = T.splitTBranches(contours)
lc = [s1]
for c in contours:
    if c[1].shape[1] == 3:
        p = G.fittingPlaster(c)
        b = G.gapfixer(c, p)
        lc.append(b)

lc = C.convertArray2Tetra(lc)
s2 = T.join(lc)
s2 = G.close(s2)
# Grille cartesienne (reguliere)
h = 0.02
ni = 200; nj = 200; nk = 2
b = G.cart( (-0.5, -0.5, 0.), (h, h, 1.), (ni,nj,nk) )
celln = C.array('cellN', ni-1, nj-1, nk-1)
celln = C.initVars(celln, 'cellN', 1.)
cellno = X.blankCells([b], [celln], [s2], blankingType=1, delta=0., dim=2)
cellno = C.center2Node(cellno, cellNType=0)
cellno[0] = T.addkplane(cellno[0])
b = C.addVars([b, cellno[0]])
# Adapte le front de la grille a la surface
b = T.subzone(b, (1,1,1), (b[2],b[3],1))
f = G.snapFront(b, [s], optimized=2)
test.testA([f], 3)
