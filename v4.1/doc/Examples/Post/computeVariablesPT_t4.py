# - computeVariables (pyTree) -
# - ro, velo, Temp
import Converter.PyTree as C
import Converter.Internal as Internal
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._initVars(a, 'centers:Density=1.')
C._initVars(a, 'centers:VelocityX=0.1')
C._initVars(a, 'centers:VelocityY=0.')
C._initVars(a, 'centers:VelocityZ=0.')
C._initVars(a, 'centers:Temperature=1.')

P._computeVariables(a, ['centers:Pressure'], rgp=(1.4-1)*1.78571)
P._computeVariables(a, ['centers:Mach'], rgp=(1.4-1)*1.78571)
P._computeVariables(a, ['centers:VelocityMagnitude'], rgp=(1.4-1)*1.78571)
P._computeVariables(a, ['centers:Enthalpy'], rgp=(1.4-1)*1.78571)
P._computeVariables(a, ['centers:Entropy'], rgp=(1.4-1)*1.78571)
P._computeVariables(a, ['centers:ViscosityMolecular'], rgp=(1.4-1)*1.78571, Ts=1, Cs=0.3831, mus=5.e-9)
P._computeVariables(a, ['centers:PressureDynamic'], rgp=(1.4-1)*1.78571)

test.testT(a, 1)
