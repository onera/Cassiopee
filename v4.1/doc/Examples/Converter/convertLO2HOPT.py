# - convertLO2HO (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D

a = D.triangle((0,0,0), (1,0,0), (1,1,0))
a = C.convertLO2HO(a, mode=0)
C.convertPyTree2File(a, 'out.cgns')
