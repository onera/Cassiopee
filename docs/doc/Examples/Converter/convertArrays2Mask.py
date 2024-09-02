# - convertArray2Mesh (Mesh Cassiopee) -
from elsA_user import *
import Converter.Cassiopee as CC
import Converter as C

# Create a cartesian mesh
msh = mesh(name='msh')
msh.submit()

# Get arrays from file
import Generator as G
a = G.cart((0.,0.,0.), (0.1,0.1,1.), (10,10,10))

# Convert arrays to mesh
CC.convertArray2Mesh(a, msh)

# Save the mesh
a = CC.convertMesh2Array(msh)
C.convertArrays2File([a], "mesh.tp")
