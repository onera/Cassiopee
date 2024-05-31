# - skySphere (pyTree) -
import Converter.PyTree as C
import Modeler.PyTree as Models

# Sky sphere
a = Models.skySphere((0,0,0), 1.)
C.convertPyTree2File(a, 'out.cgns')
