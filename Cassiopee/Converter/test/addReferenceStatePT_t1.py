# - addReferenceState (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.elsAProfile as elsAProfile
import KCore.test as test

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
t = C.newPyTree(['Base']); t[2][1][2].append(a)
tp = elsAProfile.addReferenceState(t, conservative=[1.,0,0,0,1.7])
test.testT(tp, 1)

tp = elsAProfile.addReferenceState(t, conservative=[1.,0,0,0,1.7,1e-7],turbmod='spalart')
test.testT(tp,2)
tp = elsAProfile.addReferenceState(t, conservative=[1.,0.,0.,0.,1.7,1e-7],turbmod='spalart',temp=1.,name='ReferenceStateSpalart',comments="Mon etat de reference")
test.testT(tp,3)
