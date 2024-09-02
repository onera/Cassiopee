# - computeExtraVariable (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Post.PyTree as P
import KCore.test as test

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
a = C.initVars(a, 'Density', 1.)
a = C.initVars(a, 'MomentumX', 1.)
a = C.initVars(a, 'MomentumY', 0.)
a = C.initVars(a, 'MomentumZ', 0.)
a = C.initVars(a, 'EnergyStagnationDensity', 1.)
a = P.computeExtraVariable(a, 'VorticityMagnitude')
a = P.computeExtraVariable(a, 'QCriterion')
a = P.computeExtraVariable(a, 'ShearStress')
t = C.newPyTree(['Base',a])
test.testT(t, 1)
