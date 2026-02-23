# - mergeTrees (pyTree) -
import Converter.PyTree as C

t1 = C.newPyTree(['Base1', 2, 'Base2', 3])
t2 = C.newPyTree(['Other1', 'Other2', 'Base2'])
t = C.mergeTrees(t1, t2)
C.convertPyTree2File(t, 'out.cgns')
