# - cylinder3 (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
teta = G.cart((0.,0.,0.), (0.1, 1., 1.), (11, 1, 1))
xz = G.cart((0.1,0.,0.), (0.1,1.,0.2), (20, 1, 30))
cyl = G.cylinder3(xz, 0., 90., teta)
C.convertPyTree2File(cyl, 'out.cgns')
