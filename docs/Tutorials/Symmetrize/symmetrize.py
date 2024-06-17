# Symmetrize a mesh
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.Internal as Internal

# Create a 2-block mesh for half cylinder symmetric case
a = G.cylinder((0,0,0), 0.5, 1., 180., 0., 1., (21,21,21)) 
a = C.addBC2Zone(a, 'wall', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'overlap', 'BCOverlap', 'jmax')
a = C.addBC2Zone(a, 'overlap', 'BCOverlap', 'kmin')
a = C.addBC2Zone(a, 'overlap', 'BCOverlap', 'kmax')
a1 = T.subzone(a, (1,1,1), (11,21,21))
a2 = T.subzone(a, (11,1,1), (21,21,21))
t = C.newPyTree(['MESH', a1, a2])
t = X.connectMatch(t)
t = C.fillEmptyBCWith(t, 'sym', 'BCSymmetryPlane')
t = C.initVars(t,'{F}={CoordinateX}*{CoordinateX}+{CoordinateY}+2*{CoordinateZ}')

# Symmetry y=0
t2 = T.symetrize(t, (0.,0.,0.), (1,0,0), (0,0,1))
# Rename zones with a unique name
zones = Internal.getNodesFromType(t2, 'Zone_t')
for z in zones: z[0] = C.getZoneName(z[0])
# Reorder created zones
t2 = T.reorder(t2, (-1,2,3))
# Merge both trees into one into 2 separated bases
t = C.mergeTrees(t, t2)
# Remove the BCSymmetryPlane
t = C.rmBCOfType(t, 'BCSymmetryPlane')
# Computes the BCMatch instead of symmetry plane
t = X.connectMatch(t)
# Output
C.convertPyTree2File(t, 'mesh.cgns')
