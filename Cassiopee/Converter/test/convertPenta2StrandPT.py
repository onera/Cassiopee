# - convertPenta2Strand (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cartPenta((0,0,0), (1,1,1), (3,3,3))

C._convertPenta2Strand(a)

C.convertPyTree2File(a, 'out.cgns')
