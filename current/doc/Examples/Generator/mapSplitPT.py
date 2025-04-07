# - mapSplit (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

# polyline
a = D.polyline([(0,0,0),(1,0,0),(1,1,0),(2,3,0),(1.5,3,0),(1,1.5,0),(0,0,0)])
# distribution
Ni = 41
dist = G.cart((0,0,0),(1./(Ni-1),1,1),(Ni,1,1))
dist = G.enforceX(dist, 15.5/(Ni-1), 0.005, 2,5)
dist = G.enforceX(dist, 27.5/(Ni-1), 0.005, 2,5)
zones = G.mapSplit(a,dist,0.25)
t = C.newPyTree(['Base',1]); t[2][1][2] += zones
C.convertPyTree2File(t, 'out.cgns')
