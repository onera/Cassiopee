# - computeVariables (array) -
import Converter as C
import Post as P
import Generator as G
import KCore.test as test

# test sur un array
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
c = C.array('ro,rou,rov,row,roE', ni, nj, 2)
c = C.initVars(c, 'ro', 1.)
c = C.initVars(c, 'rou', 1.)
c = C.initVars(c, 'rov', 0.)
c = C.initVars(c, 'row', 0.)
c = C.initVars(c, 'roE', 1.)
m = C.addVars([m,c])
A = [m]
vars = ['Pressure', 'Mach', 'Entropy', 'Enthalpy',
        'VelocityX', 'VelocityY', 'VelocityZ',
        'Temperature', 'ViscosityMolecular',
        'PressureStagnation', 'TemperatureStagnation',
        'PressureDynamic']
sol = P.computeVariables(m, vars)
test.testA([sol], 1)

# test sur une liste
ni = 20; nj = 10
m = G.cart((0,0,0), (5./(ni-1),10./(nj-1),1), (ni,nj,2))
c = C.array('ro,rou, rov,row,roE', ni, nj, 2)
c = C.initVars(c, 'ro', 1.)
c = C.initVars(c, 'rou', 1.)
c = C.initVars(c, 'rov', 0.)
c = C.initVars(c, 'row', 0.)
c = C.initVars(c, 'roE', 1.)
m = C.addVars([m,c])
A.append(m)
vars = ['Pressure', 'Mach']
sol = P.computeVariables(A, vars)
test.testA(sol, 2)
#
vars = ['PressureStagnation','TemperatureStagnation','PressureDynamic']
sol = P.computeVariables(A, vars)
test.testA(sol, 3)

# Test avec des grandeurs ro,u,T
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
c = C.array('ro,u,v,w,T', ni, nj, 2)
c = C.initVars(c, 'ro', 1.)
c = C.initVars(c, 'u', 1.)
c = C.initVars(c, 'v', 0.)
c = C.initVars(c, 'w', 0.)
c = C.initVars(c, 'Temp', 1.)
m = C.addVars([m,c])
A = [m]
vars = ['Pressure', 'Mach', 'Entropy', 'Enthalpy',
        'VelocityX', 'VelocityY', 'VelocityZ',
        'Temperature', 'ViscosityMolecular',
        'PressureStagnation', 'TemperatureStagnation',
        'PressureDynamic']
sol = P.computeVariables(m, vars)
test.testA([sol], 4)
