# Computes the matching connectivity for a portion of cylinder
import Converter.PyTree as C
import Connector.PyTree as X
import Transform.PyTree as T
import Generator.PyTree as G

a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 90., 5., (11,11,11)) 
t = C.newPyTree(['Cylindre', a])

# Duplicate mesh by a rotation of +90 degrees
ap = T.rotate(a, (0.,0.,0.),(0.,0.,1.), 90.)
t = C.addBase2PyTree(t, 'CylindreDup'); t[2][2][2] += [ap]
t = X.connectMatch(t)

# Duplicate mesh by a rotation of +90 degrees
# Replace duplicated zones
am = T.rotate(a, (0.,0.,0.), (0.,0.,1.), -90.) ; t[2][2][2] = [am]
t = X.connectMatch(t)

# Remove duplicated basis
t[2][1:] = t[2][1:2]
C.convertPyTree2File(t, 'out.cgns')
