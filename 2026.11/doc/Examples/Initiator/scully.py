# - initScully (array) -
import Generator as G
import Converter as C
import Initiator as I

NI = 200; NJ = 200
HI = 1./(NI-1); HJ = 1./(NJ-1)
a = G.cart((0.,0.,0.), (HI,HJ,1.), (NI,NJ,2))
a = I.initScully(a, (0.5,0.5), -0.2, 0.05, 0.8, 0)
C.convertArrays2File([a], "out.plt")
