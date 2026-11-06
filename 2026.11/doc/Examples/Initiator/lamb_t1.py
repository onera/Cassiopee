# - initLamb -
import Generator as G
import Initiator as I
import KCore.test as test

# Structure
NI = 100; NJ = 100
HI = 50./(NI-1); HJ = 50./(NJ-1)
a = G.cart( (0.,0.,0.), (HI,HJ,1.), (NI,NJ,2))
a = I.initLamb(a, position=(25.,25.), Gamma=2., MInf=0.8)
test.testA([a], 1)

# Non structure
NI = 100; NJ = 100
HI = 50./(NI-1); HJ = 50./(NJ-1)
a = G.cartTetra( (0.,0.,0.), (HI,HJ,1.), (NI,NJ,2))
a = I.initLamb(a, position=(25.,25.), Gamma=2., MInf=0.8)
test.testA([a], 2)
