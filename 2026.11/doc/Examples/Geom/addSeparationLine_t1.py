# - addSeparationLine (array) -
import Geom as D
import Converter as C
import KCore.test as test

# Add a line to a circle
a1 = D.circle((0,0,0), 1, 0., 360, 1000)
a2 = D.line((0.,1.,0.), (0.,2.,0), 100)
a0 = D.addSeparationLine(a1, a2)
test.testA(a0, 1)

# Avec un demi cercle
a1 = D.circle((0,0,0), 1, 180, 360, 1000)
a2 = D.line((0.,-1.,0.), (0.,-2.,0), 100)
a0 = D.addSeparationLine(a1, a2)
test.testA(a0, 2)

# Avec un champ en plus
a1 = D.circle((0,0,0), 1, 0., 360, 1000)
a1 = C.initVars(a1, '{F}=3*{x}+{y}')
a2 = D.line((0.,1.,0.), (0.,2.,0), 100)
a2 = C.initVars(a2, '{F}=3*{x}+{y}')
a0 = D.addSeparationLine(a1, a2)
test.testA(a0, 3)
