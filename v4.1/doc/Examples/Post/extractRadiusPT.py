import Transform.PyTree as T
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Post.Rotor as Rotor
import math

a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,30))
coordZ = Internal.getNodeFromName(a, 'CoordinateZ')[1]
coordX = Internal.getNodeFromName(a, 'CoordinateX')[1]
Internal.getNodeFromName(a, 'CoordinateZ')[1] = coordX
Internal.getNodeFromName(a, 'CoordinateX')[1] = coordZ

axis_pnt = (0,0,0)
axis_vct = [1,0,0]

theta = 25.
costheta = math.cos(math.radians(theta))
sintheta = - math.sin(math.radians(theta))
a = T.rotate(a, (0,0,0), (0.,theta,0.))
a = T.rotate(a, (0,0,0), (0.,0.,theta))
a = T.rotate(a, (0,0,0), (theta,0.,0.))

#rotation y
x,y,z = axis_vct
axis_vct[0] = costheta*x + sintheta*z
axis_vct[1] = y
axis_vct[2] = - sintheta*x + costheta*z

#rotation z
x,y,z = axis_vct
axis_vct[0] = costheta*x + sintheta*y
axis_vct[1] = - sintheta*x + costheta*y
axis_vct[2] = z

#rotation x
x,y,z = axis_vct
axis_vct[0] = x
axis_vct[1] = costheta*y + sintheta*z
axis_vct[2] = - sintheta*y + costheta*z

a = Rotor.extractRadius(a, axis_pnt, axis_vct)

C.convertPyTree2File(a, 'out.cgns')