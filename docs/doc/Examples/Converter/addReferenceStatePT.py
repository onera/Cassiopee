# - addReferenceState (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.elsAProfile as elsAProfile

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a])
tp = elsAProfile.addReferenceState(t, conservative=[1.,0,0,0,1.7])
C.convertPyTree2File(tp, 'out.cgns')
