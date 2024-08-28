# - initLamb (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Initiator.PyTree as I

NI = 100; NJ = 100
HI = 50./(NI-1); HJ = 50./(NJ-1)
a = G.cart((0.,0.,0.), (HI,HJ,1.), (NI,NJ,2))
a = I.initLamb(a, position=(7.,7.), Gamma=2., MInf=0.8, loc='centers')
C.convertPyTree2File(a, "out.cgns")
