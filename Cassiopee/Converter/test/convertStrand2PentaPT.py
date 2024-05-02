# - convertStrand2Penta (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cartPenta((0,0,0), (1,1,1), (3,3,3))

b = C.convertPenta2Strand(a)

C._convertStrand2Penta(b)

C.convertPyTree2File(b, 'out.cgns')
