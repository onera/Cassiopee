# - initConst (array) -
import KCore.test as test
import Generator as G
import Initiator as I

# Structure
NI = 100; NJ = 100
HI = 50./(NI-1); HJ = 50./(NJ-1)
a = G.cart( (0.,0.,0.), (HI,HJ,1.), (NI,NJ,2))
a = I.initConst(a, MInf=0.8)
test.testA([a], 1)

# Non structure
NI = 100; NJ = 100
HI = 50./(NI-1); HJ = 50./(NJ-1)
a = G.cartTetra( (0.,0.,0.), (HI,HJ,1.), (NI,NJ,2))
a = I.initConst(a, MInf=0.8)
test.testA([a], 2)
