# - getCurvilinearAbscissa (pyTree)-
import Converter.PyTree as C
import Geom.PyTree as D

a = D.line((0.,0.,0.), (1.,0.,0), 100)
a = D.getCurvilinearAbscissa(a)
C.convertPyTree2File(a, 'out.cgns')
