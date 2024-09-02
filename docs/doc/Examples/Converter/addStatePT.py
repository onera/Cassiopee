# - addState (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cylinder((0,0,0), 1., 1.5, 0., 360., 1., (80,30,2))
t = C.newPyTree(['Base',a])

# Specifie des valeurs
b = Internal.getNodeFromName1(t, 'Base')
C._addState(b, 'EquationDimension', 2)
C._addState(b, 'GoverningEquations', 'Euler')
C._addState(b, 'Mach', 0.6)
C._addState(b, 'Reynolds', 100000)
