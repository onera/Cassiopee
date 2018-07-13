# - overlayField (pyTree) -
import Generator.PyTree as G
import Initiator.PyTree as I
import Converter.PyTree as C
import KCore.test as test

NI = 100; NJ = 100
HI = 50./(NI-1); HJ = 50./(NJ-1)
MachInf = 0.8
a = G.cart( (0.,0.,0.), (HI,HJ,1.), (NI,NJ,2))
a = C.addBC2Zone(a, 'wall','BCWall','imin')
a = C.addBC2Zone(a, 'overlap', 'BCOverlap', 'jmin')
z1 = I.initVisbal(a, (3.,3.), 2., MachInf, loc='centers'); z1[0]='cart1'
z2 = I.initVisbal(a, (20.,3.), 2., MachInf, loc='centers'); z2[0]='cart2'
an = I.overlayField(z1, z2, MachInf, loc='centers')
t = C.newPyTree(['Base']); t[2][1][2] += [an]
test.testT(t,1)
