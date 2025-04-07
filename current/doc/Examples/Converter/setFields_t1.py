# - setField (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
a = C.initVars(a, "F", 1.); a = C.initVars(a, "centers:Q", 1.2)
t = C.newPyTree(['Cart',a])
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
ac = C.getFields('FlowSolution#Centers', t)
t = C.setFields(ac, t, 'centers')
test.testT(t)
