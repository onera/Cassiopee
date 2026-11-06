# - computeExtraVariable (array) -
import Generator as G
import Converter as C
import Post as P
import Transform as T
import KCore.test as test

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
a = C.initVars(a, 'Density', 1.)
a = C.initVars(a, 'MomentumX', 1.)
a = C.initVars(a, 'MomentumY', 0.)
a = C.initVars(a, 'MomentumZ', 0.)
a = C.initVars(a, 'EnergyStagnationDensity', 1.)
v = P.computeExtraVariable(a, 'Vorticity')
m = P.computeExtraVariable(a, 'VorticityMagnitude')
q = P.computeExtraVariable(a, 'QCriterion')
tau = P.computeExtraVariable(a, 'ShearStress')
a = C.node2Center(a)
a = C.addVars([a, v, m, q, tau])
test.testA([a], 1)
wall = T.subzone(a, (1,1,1), (9,9,1))
skinFriction = P.computeExtraVariable(wall, 'SkinFriction')
wall = C.addVars([wall, skinFriction])
test.testA([wall], 2)
