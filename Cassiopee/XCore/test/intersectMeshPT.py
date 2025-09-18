# - intersectMesh (pyTree) -
import Converter.PyTree as C
import Converter.Internal as I
import XCore.PyTree as X
import Generator.PyTree as G
import Transform.PyTree as T

# Load master hexa mesh and transform it to NGon
m = G.cartNGon((0,0,0),(0.1,0.1,0.1),(11,11,11))
T._rotate(m,(0,0,0),(5.0,1.0,2.))
# Mark the original cells to keep
C._initVars(m, 'centers:keep', 1.0)
# Triangulate the external quads
X._triangulateSkin(m)

# Initialize the IntersectMesh hook
IM = X.IntersectMesh_Init(m)

# Load the structured meshes to modify
c = G.cart((0.19,0.17,-0.2),(0.05,0.05,0.05),(15,15,15))
c = C.newPyTree(['cylinder', c])
c = G.close(c)
T._rotate(c,(0,0,0),(10,2.0,4.))

# Remove the K-planes within the master mesh, and project
cp = X.removeIntersectingKPlanes(IM, c)
# Mark the original cells to keep
C._initVars(cp, 'centers:keep', 1.0)

# Drop the IntersectMesh hook
X.IntersectMesh_Exit(IM)

# Initialize the adaptation-intersection hook
IC = X.icapsuleInit2()
# Set the master mesh
X.icapsuleSetMaster(IC, m)
# Set the slave meshes
S = C.newPyTree(["Slave"]); S[2][1][2]=I.getZones(cp)
X.icapsuleSetSlaves(IC, [S])

# Intersect
X.icapsuleAdapt2(IC)
X.icapsuleIntersect2(IC)
# Extract the resulting meshes as zones
Mi = X.icapsuleExtractMaster(IC)
Si = X.icapsuleExtractSlaves(IC)
C.convertPyTree2File(Mi, "Mi.cgns")
C.convertPyTree2File(Si, "Si.cgns")
