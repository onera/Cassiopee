# - overlayField (array) -
import Converter as C
import Generator as G
import Initiator as I

NI = 100; NJ = 100
HI = 50./(NI-1); HJ = 50./(NJ-1)
MachInf = 0.8

a = G.cart((0.,0.,0.), (HI,HJ,1.), (NI,NJ,2))
ac = I.initVisbal(a, (3.,3.), 2., MachInf)
ac2 = I.initVisbal(a, (20.,3.), 2., MachInf)
an = I.overlayField(ac, ac2, MachInf)
C.convertArrays2File([an], "out.plt")
