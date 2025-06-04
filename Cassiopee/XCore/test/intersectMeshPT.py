# - intersectMesh (pyTree) -
import Converter.PyTree as C
import Converter.Internal as I
import XCore.PyTree as X
import Generator.PyTree as G

# Load master hexa mesh and transform it to NGon
m = G.cartNGon((0,0,0),(0.1,0.1,0.1),(11,11,11))
# Mark the original cells to keep
C._initVars(m, 'centers:keep', 1.0)
# Triangulate the external quads
X._triangulateSkin(m)
I._adaptNGon42NGon3(m)

# Initialize the IntersectMesh hook
IM = X.IntersectMesh_Init(m)

# Load the structured meshes to modify
c = G.cart((0.19,0.17,-0.5),(0.11,0.12,0.1),(5,5,5))
c = C.newPyTree(['cylinder', c])
c = G.close(c)

# Remove the K-planes within the master mesh, and project
cp = X.removeIntersectingKPlanes(IM, c)
# Mark the original cells to keep
C._initVars(cp, 'centers:keep', 1.0)
I._adaptNGon42NGon3(cp)

# Drop the IntersectMesh hook
X.IntersectMesh_Exit(IM)

# Initialize the adaptation-intersection hook
IC = X.icapsuleInit2()
# Set the master mesh
X.icapsuleSetMaster(IC, m)
# Set the slave meshes
X.icapsuleSetSlaves(IC, [cp])

# Intersect
X.icapsuleIntersect2(IC)
# Extract the resulting meshes as zones
Mi = X.icapsuleExtractMaster(IC)
Si = X.icapsuleExtractSlaves(IC)
C.convertPyTree2File(Mi, "Mi.cgns")
C.convertPyTree2File(Si, "Si.cgns")
