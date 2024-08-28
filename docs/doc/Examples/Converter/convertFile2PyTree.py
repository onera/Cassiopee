# - convertFile2PyTree (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
t = C.newPyTree(['Base',a])
C.convertPyTree2File(t, 'in.cgns')
t1 = C.convertFile2PyTree('in.cgns'); print(t1)
