# - addNormalLayers (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

d = G.cart((0.1,0.,0.), (0.1,1,1),(2,1,1)) 
a = D.sphere((0,0,0), 1, 50)
a = D.line((0,0,0),(1,1,0),3)
b = G.addNormalLayers(a, d)
d = G.cart((0.1,0.,0.), (-0.1,1,1),(2,1,1)) 
c =  G.addNormalLayers(a, d)
import Transform.PyTree as T
a = T.join(b,c)
C.convertPyTree2File(a, 'out.cgns')
