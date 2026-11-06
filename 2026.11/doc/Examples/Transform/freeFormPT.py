# - freeForm (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import Transform.PyTree as T

# shape
a = D.circle((0,0,0), 1., N=50)

# create control points
b = T.controlPoints(a, (2,2,1))
C.setValue(b, 'dx', (0,0,0), -0.2)
C.setValue(b, 'dy', (0,0,0), -0.2)

# free form
T._freeForm(a, b)

a = T.deform(a, ['dx','dy','dz'])
b = T.deform(b, ['dx','dy','dz'])

C.convertPyTree2File([a,b], 'out.cgns')
