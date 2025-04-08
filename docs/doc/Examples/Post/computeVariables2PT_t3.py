# - computeVariables2PT (pyTree) -
import Converter.PyTree  as CP
import Post.PyTree       as PT
import Generator.PyTree  as GP
import KCore.test        as test

ni = 30; nj = 40
z1 = GP.cartNGon((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))

# Variables a calculer
# --------------------
vars = ['centers:Pressure','VelocityX','VelocityZ','VelocityMagnitude',
        'Temperature','Entropy','centers:Enthalpy','Mach','ViscosityMolecular',
        'PressureStagnation','centers:TemperatureStagnation']

# Variables conservatives
# -----------------------
z1 = CP.initVars(z1, 'Density', 1.)
z1 = CP.initVars(z1, 'MomentumX', 1.)
z1 = CP.initVars(z1, 'MomentumY', 0.)
z1 = CP.initVars(z1, 'MomentumZ', 0.)
z1 = CP.initVars(z1, 'EnergyStagnationDensity', 1.e5)
z1 = CP.initVars(z1, 'Mach', 3.)

z1 = CP.initVars(z1, 'centers:Density', 5.)
z1 = CP.initVars(z1, 'centers:MomentumX', 15.)
z1 = CP.initVars(z1, 'centers:MomentumY', 3.)
z1 = CP.initVars(z1, 'centers:MomentumZ', 2.)
z1 = CP.initVars(z1, 'centers:EnergyStagnationDensity', 2.e5)

# computeVariables2 sur une zone
# ------------------------------
z1 = PT.computeVariables2(z1, vars)
test.testT(z1,1)

# computeVariables2 sur un PyTree
# -------------------------------
t1 = CP.newPyTree(['Base', z1])
t1 = CP.addState(t1, 'Gamma', 5.9)
t1 = PT.computeVariables2(t1, vars)

test.testT(t1,2)

# # Variables primitives
# # ---------------------
z2 = GP.cartNGon((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))

z2 = CP.initVars(z2, 'Density', 1.)
z2 = CP.initVars(z2, 'VelocityX', 1.)
z2 = CP.initVars(z2, 'VelocityY', 0.)
z2 = CP.initVars(z2, 'VelocityZ', 0.)
z2 = CP.initVars(z2, 'Temperature', 1.)

z2 = CP.initVars(z2, 'centers:Density', 1.)
z2 = CP.initVars(z2, 'centers:VelocityX', 1.)
z2 = CP.initVars(z2, 'centers:VelocityY', 0.)
z2 = CP.initVars(z2, 'centers:VelocityZ', 0.)
z2 = CP.initVars(z2, 'centers:Temperature', 1.)

# computeVariables2 sur une zone
# ------------------------------
PT._computeVariables2(z2, vars)
test.testT(z2,3)

# # computeVariables2 sur un PyTree
# # -------------------------------
t2 = CP.newPyTree(['Base', z2])
t2 = CP.addState(t2, 'Cv', 10.9/0.4)
PT._computeVariables2(t2, vars)
test.testT(t2,4)
