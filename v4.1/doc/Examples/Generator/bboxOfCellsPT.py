# - bboxOfCells (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(20,20,20))
a = G.bboxOfCells(a)
C.convertPyTree2File(a, 'out.cgns')
