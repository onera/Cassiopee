# - setField (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
a = C.initVars(a, "F", 1.); a = C.initVars(a, "centers:Q", 1.2)
t = C.newPyTree(['Cart',a])
ac = C.getFields('FlowSolution#Centers', t)
t = C.setFields(ac, t, 'centers')
t = C.addBase2PyTree(t, 'Extra'); t[2][2][2].append(a)
C.convertPyTree2File(t, 'out.cgns')
