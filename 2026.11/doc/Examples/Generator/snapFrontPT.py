# - snapFront (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Connector.PyTree as X
import Transform.PyTree as T
import Converter.Internal as Internal

s = D.circle((0,0,0), 1., N=100)
s2 = T.addkplane(s)

# Grille cartesienne (reguliere)
BB = G.bbox([s])
ni = 100; nj = 100; nk = 3
xmin = BB[0]; ymin = BB[1]; zmin = BB[2]-0.5
xmax = BB[3]; ymax = BB[4]; zmax = BB[5]+0.5
hi = (xmax-xmin)/(ni-1); hj = (ymax-ymin)/(nj-1)
h = min(hi, hj)
ni = int((xmax-xmin)/h)+7; nj = int((ymax-ymin)/h)+7
b = G.cart( (xmin-3*h, ymin-3*h, zmin), (h, h, 1.), (ni,nj,nk) )
t = C.newPyTree(['Cart'])
t[2][1][2].append(b)

# Masquage
t = C.initVars(t,'cellN',1)
import numpy
BM = numpy.array([[1]])
t = X.blankCells(t,[[s2]],BM,blankingType='node_in',dim=2)

# Adapte le front de la grille a la surface
dim = Internal.getZoneDim(b)
t = T.subzone(t, (1,1,2), (dim[1],dim[2],2))
t = G.snapFront(t, [s2])

t = C.addBase2PyTree(t, 'Surface', cellDim=2)
s2 = C.initVars(s2,'cellN',1)
t[2][2][2].append(s2)

C.convertPyTree2File(t, 'out.cgns')
