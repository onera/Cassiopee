# - cone (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import KCore.test as test

a = D.cone((0,0,0), 1. , 0.5, 1.)
t = C.newPyTree(['Base',2,a])
test.testT(t, 1)
