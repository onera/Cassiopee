# - refine (array) -
import Post as P
import Converter as C
import Generator as G
import KCore.test as test

a = G.cartTetra( (0,0,0), (1,1,1), (10,10,1))
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
indic = C.array('indic', a[2].shape[1], 1, 1)
indic = C.initVars(indic, 'indic', 0)
C.setValue(indic, 50, [1])
C.setValue(indic, 49, [1])
a = P.refine(a, indic)
test.testA([a])
