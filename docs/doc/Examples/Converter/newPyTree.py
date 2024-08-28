# - newPyTree (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal

# Create a tree with two bases
t = C.newPyTree(['Base1','Base2'])

# Create a tree with two bases with their dims
t = C.newPyTree(['Base1',2,'Base2',3])

# Create a tree with an attached existing Base node
base = Internal.newCGNSBase('Base', 3)
t = C.newPyTree([base])

# Create a tree with existing zones attached
z1 = Internal.newZone('Zone1')
z2 = Internal.newZone('Zone2')
t1 = C.newPyTree(['Base', z1, z2])
t2 = C.newPyTree(['Base', z1, 'Base2', z2])
t3 = C.newPyTree(['Base', [z1,z2]])
C.convertPyTree2File(t3, 'out.cgns')
