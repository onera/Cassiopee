# - map (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.PyTree as C

l = D.line( (0,0,0), (1,1,0) )
Ni = 11; dist = G.cart( (0,0,0), (1./(Ni-1),1.,1.), (Ni,1,1) )
l = G.map(l, dist)
t = C.newPyTree(['Base',1,l])
C.convertPyTree2File(t, 'out.cgns')
