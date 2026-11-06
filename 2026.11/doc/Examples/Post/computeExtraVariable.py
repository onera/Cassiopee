# - computeExtraVariable (array) -
import Generator as G
import Converter as C
import Post as P
import Transform as T

a = G.cart( (0,0,0), (1,1,1), (50,50,50) )
a = C.initVars(a, 'Density', 1.)
a = C.initVars(a, 'MomentumX', 1.)
a = C.initVars(a, 'MomentumY', 0.)
a = C.initVars(a, 'MomentumZ', 0.)
a = C.initVars(a, 'EnergyStagnationDensity', 100000.)
m = P.computeExtraVariable(a, 'VorticityMagnitude')
q = P.computeExtraVariable(a, 'QCriterion')
tau = P.computeExtraVariable(a, 'ShearStress')
a = C.node2Center(a)
a = C.addVars([a, m, q, tau])

# Skin friction requires a surface array with shear stress already computed
wall = T.subzone(a, (1,1,1), (49,49,1))
skinFriction = P.computeExtraVariable(wall, 'SkinFriction')
skinFrictionTangential = P.computeExtraVariable(wall, 'SkinFrictionTangential')
wall = C.addVars([wall, skinFriction, skinFrictionTangential])
C.convertArrays2File([wall], 'out.plt')
