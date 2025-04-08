# - initVisbal -
import Generator as G
import Initiator as I
import KCore.test as test
import Converter as C

NI = 100; NJ = 100
HI = 50./(NI-1); HJ = 50./(NJ-1)
MachInf = 0.8

# Structure
a = G.cart( (0.,0.,0.), (HI,HJ,1.), (NI,NJ,2))
a = C.initVars(a, 'Density', 1.)
ac = I.initVisbal(a, (25.,25.), 2., 0.8)
test.testA([ac],1)

# Non structure
a = G.cartTetra( (0.,0.,0.), (HI,HJ,1.), (NI,NJ,2))
a = C.initVars(a, 'Density', 1.)
ac = I.initVisbal(a, (25.,25.), 2., 0.8)
test.testA([ac],2)
