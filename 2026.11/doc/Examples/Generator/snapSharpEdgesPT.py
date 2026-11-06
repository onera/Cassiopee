# - snapSharpEdges (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

# polylignes avec angles vifs
s = D.polyline([(0.02,0,0),(1,1,0),(2,1,0),(0.02,0,0)])
# Grille cartesienne (reguliere)
h = 0.1
ni = 30; nj = 20; nk=1
b = G.cart((-0.5, -0.5, 0), (h, h, 1.), (ni,nj,nk))
b = G.snapSharpEdges(b, [s], h)
t = C.newPyTree(['Cart','Surface'])
t[2][1][2].append(b); t[2][2][2].append(s)
C.convertPyTree2File(t, 'out.cgns')
