# - torus (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.torus((0.,0.,0.), 5., 2.)
C.convertPyTree2File(a, 'out.cgns')
