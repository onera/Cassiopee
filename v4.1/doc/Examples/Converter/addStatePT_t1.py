# - addState (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cylinder((0,0,0), 1., 1.5, 0., 360., 1., (80,30,2))
a = C.addBC2Zone(a,'wall','BCWall','imin')
a = C.addBC2Zone(a,'overlap','BCOverlap','imax')
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)

t = C.newPyTree(['Base',3]); t[2][1][2].append(a)
# Add state by value
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t[2][1] = C.addState(t[2][1], 'GoverningEquations', 'Euler')
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t[2][1] = C.addState(t[2][1], 'Reynolds', 100000.)
test.testT(t, 1)

# Add by ref state adim1
t = C.newPyTree(['Base',3]); t[2][1][2].append(a)
t = C.addState(t, adim='adim1', MInf=0.5)
test.testT(t, 2)

# Add by ref state adim2
t = C.newPyTree(['Base',3]); t[2][1][2].append(a)
t = C.addState(t, adim='adim2', MInf=0.5)
test.testT(t, 3)

# Add by ref state dim2
t = C.newPyTree(['Base',3]); t[2][1][2].append(a)
t = C.addState(t, adim='dim2', UInf=35., TInf=294., RoInf=1.2, LInf=1.)
test.testT(t, 3)
