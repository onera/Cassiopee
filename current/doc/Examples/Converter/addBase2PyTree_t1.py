# - addBase2PyTree (pyTree) -
import Converter.PyTree as C
import KCore.test as test

# Sur un arbre existant
t = C.newPyTree(['Base', 3])
t = C.addBase2PyTree(t, 'Base2', 2)
C._addBase2PyTree(t, 'Base3', 3)
test.testT(t,1)

# Sur un arbre vide
t2 = C.addBase2PyTree([], 'Base3', 3)
test.testT(t2,2)

test.writeCoverage(100.)
