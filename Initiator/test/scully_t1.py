# - initScully (array)-
import Generator as G
import Initiator as I
import KCore.test as test

# Structure
NI = 100; NJ = 100
HI = 1./(NI-1); HJ = 1./(NJ-1)
a = G.cart((0.,0.,0.), (HI,HJ,1.), (NI,NJ,2))
# constant entropy outward + constant temperature in the vortex
ac = I.initScully(a, (0.5,0.5), -0.2, 0.05, 0.8, 0)
test.testA([ac], 1)
# constant entropy
ac = I.initScully(a, (0.5,0.5), -0.2, 0.05, 0.8, 1)
test.testA([ac], 2)

# Non structure
NI = 100; NJ = 100
HI = 1./(NI-1); HJ = 1./(NJ-1)
a = G.cartTetra((0.,0.,0.), (HI,HJ,1.), (NI,NJ,2))
# constant entropy outward + constant temperature in the vortex
ac = I.initScully(a, (0.5,0.5), -0.2, 0.05, 0.8, 0)
test.testA([ac], 3)
# constant entropy
ac = I.initScully(a, (0.5,0.5), -0.2, 0.05, 0.8, 1)
test.testA([ac], 4)
