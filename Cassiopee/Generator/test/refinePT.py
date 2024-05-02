# - refine (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cart( (0,0,0), (0.1,0.1,0.1), (20,20,1) )
a = G.refine(a, 1.5, 1)
C.convertPyTree2File(a, 'out.cgns')
