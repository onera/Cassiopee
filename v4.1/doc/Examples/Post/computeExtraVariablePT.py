# - computeExtraVariable (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Transform.PyTree as T
import Post.PyTree as P

def F(x,y): return x*x + y*y

a = G.cart((0,0,0), (1,1,1), (50,50,50))
a = C.initVars(a, 'Density', 1.)
a = C.initVars(a, 'MomentumX', F, ['CoordinateX', 'CoordinateY'])
a = C.initVars(a, 'MomentumY', 0.)
a = C.initVars(a, 'MomentumZ', 0.)
a = C.initVars(a, 'EnergyStagnationDensity', 100000.)
a = P.computeExtraVariable(a, 'centers:VorticityMagnitude')
a = P.computeExtraVariable(a, 'centers:QCriterion')
a = P.computeExtraVariable(a, 'centers:ShearStress')

b = T.subzone(a, (1,1,1), (50,50,1))
b = P.computeExtraVariable(b, 'centers:SkinFriction')
b = P.computeExtraVariable(b, 'centers:SkinFrictionTangential')

C.convertPyTree2File(a, 'out.cgns')
