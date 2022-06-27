# - extractPressure (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Initiator.PyTree as I
import Post.ExtraVariables2 as PE

a = G.cart((0.,0.,0.), (13./100.,13./100.,1.), (100,100,2))
I._initLamb(a, position=(7.,7.), Gamma=2., MInf=0.8, loc='centers')
I._cons2Prim(a)
PE._extractPressure(a)
C.convertPyTree2File(a, 'out.cgns')
