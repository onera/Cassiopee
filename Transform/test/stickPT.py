# - stick (pyTree) -
import Geom.PyTree as D
import Transform.PyTree as T
import Converter.PyTree as C
import Generator.PyTree as G

# Pale
naca = D.naca(12.)
naca = T.reorder(naca, (-1,2,3))
line = D.line((0,0,0), (0,0,1))
a = D.lineDrive(naca, line)
a = T.rotate(a, (0,0,0), (1,0,0), 90.)
d = D.line((0,0,0),(0.1,0,0),3)
a = G.addNormalLayers(a, d)
C._addBC2Zone(a, 'stick', 'FamilySpecified:stick', 'jmin')
C._fillEmptyBCWith(a, 'free', 'FamilySpecified:free')
t = C.newPyTree(['Base', a])
C.convertPyTree2File(t, 'pale.cgns')

# Hub
circle = D.circle( (0.5,0.7,0), 0.9 )
line = D.line((0,0,0), (0,0,0.5))
hub = D.lineDrive(circle, line)
hub = T.translate(hub, (0,0,-0.25))
C.convertPyTree2File(hub, 'hub.cgns')

# stick
T._stick(t, hub)
C.convertPyTree2File(t, 'out.cgns')

