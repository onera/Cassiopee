# - initYee (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Initiator.PyTree as I
import KCore.test as test

NI = 100; NJ = 100
HI = 50./(NI-1); HJ = 50./(NJ-1)
MachInf = 0.8

a = G.cart( (0.,0.,0.), (HI,HJ,1.), (NI,NJ,2))
z = I.initYee(a, (3.,3.), 2., MachInf, loc='centers')
t = C.newPyTree(['Base']); t[2][1][2].append(z)
test.testT(t, 1)
