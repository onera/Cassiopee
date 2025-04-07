import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

# Create a cartesian grid
a = G.cart( (0,0,0), (1,1,1), (10,11,12) )

# A node is a list of type: ['Name', n, [], 'NodeType_t']
# 'Name' is the name if the node, n is a numpy containing the data,
# [] is a list of sons of this node, 'NodeType_t' describes the type of node.

# This is a zone node
Internal.printTree(a)

# To manipulate nodes, we use the Internal module
node = Internal.getNodeFromName(a, 'CoordinateX')

# n is a numpy you can manipulate classicaly
n = node[1]

# Convert a zone in unstructured HEXA
# This function returns a copy of the zone
b = C.convertArray2Hexa(a)

# You can directly modify a without a copy (in place functions), using
# the function prefixed with a _
C._convertArray2Hexa(a)

# You can save as a pyTree
C.convertPyTree2File([a], 'out.cgns')
