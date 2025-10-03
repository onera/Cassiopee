# - freeForm (array) -
import Converter as C
import Geom as D
import Transform as T

# shape
a = D.circle((0,0,0), 1., N=50)

# create control points
b = T.controlPoints(a, (2,2,1))
b[1][3,0] = -0.2; b[1][4,0] = -0.2 # deformation
#b[1][3,:] = 0.2 # translation globale

# free form
a = T.freeForm(a, b)

a = T.deform(a, ['dx','dy','dz'])
b = T.deform(b, ['dx','dy','dz'])

C.convertArrays2File([a,b], 'out.plt')
