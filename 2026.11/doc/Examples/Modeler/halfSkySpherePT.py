# - halfSkySphere (pyTree) -
import Converter.PyTree as C
import Modeler.PyTree as Models

# Half sky sphere
a = Models.halfSkySphere((0,0,0), 1.)
C.convertPyTree2File(a, 'out.cgns')
