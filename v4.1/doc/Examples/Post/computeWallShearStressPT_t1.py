# - computeWallShearStress (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Post.PyTree as P
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (50,50,1))
t = C.newPyTree(['Base',a])
C._addState(t, state='EquationDimension', value=3)
C._addState(t, adim='adim1')
C._initVars(t,'{VelocityX}=0.2*{CoordinateX}**2')
C._initVars(t,'{VelocityY}=0.3*{CoordinateY}*{CoordinateX}')
C._initVars(t,'VelocityZ', 0.)
for var in ['VelocityX','VelocityY','VelocityZ']:
    t = P.computeGrad(t,var)
    t = C.node2Center(t,var)
C._initVars(t,'centers:Density', 1.)
C._initVars(t,'centers:EnergyStagnationDensity', 1.)
C._initVars(t,'centers:Temperature', 1.)
tw = P.computeWallShearStress(t)
test.testT(tw,1)
#
t = C.newPyTree(['Base',a])
C._addState(t, state='EquationDimension', value=3)
C._addState(t, adim='adim1')
C._initVars(t,'{VelocityX}=0.2*{CoordinateX}**2')
C._initVars(t,'{VelocityY}=0.3*{CoordinateY}*{CoordinateX}')
C._initVars(t,'VelocityZ', 0.)
for var in ['VelocityX','VelocityY','VelocityZ']:
    t = P.computeGrad(t,var)
    t = C.center2Node(t,['centers:gradx'+var,'centers:grady'+var,'centers:gradz'+var])
C._initVars(t,'Density', 1.)
C._initVars(t,'EnergyStagnationDensity', 1.)
C._initVars(t,'Temperature', 1.)
tw = P.computeWallShearStress(t)
test.testT(tw,2)
