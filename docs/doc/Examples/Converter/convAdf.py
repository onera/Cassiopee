# - convAdf (pyTree) -
import Converter.PyTree as C

t = C.convertFile2PyTree('out.cgns')
C.convertPyTree2File(t, 'out2.cgns')
