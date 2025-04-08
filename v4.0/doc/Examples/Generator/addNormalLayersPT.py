# - addNormalLayers (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

d = G.cart((0.,0.,0.), (0.1,1,1),(3,1,1))
a = D.sphere((0,0,0), 1, 50)
b = G.addNormalLayers(a, d)
C.convertPyTree2File(b, 'out.cgns')
