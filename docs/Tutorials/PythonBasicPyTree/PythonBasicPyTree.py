import Converter.PyTree as C
import Converter.Internal as Internal

# Read file
a = C.convertFile2PyTree("case.cgns")
# print pytree a to screen
Internal.printTree(a)

# get a node by its name (first met)
n1 = Internal.getNodeFromName(a, 'CoordinateX')
# get a node by its type (first met)
n2 = Internal.getNodeFromType(a, 'FlowSolution_t')
# get all nodes by its type (first met)
n3 = Internal.getNodesFromType(a, 'FlowSolution_t')
# you can start search from any node
n4 = Internal.getNodeFromName(n3, 'Density')

# get numpy array of tree node
array = n1[1]
# change the X coordinate of first point
array[0] = 0.1

# save pytree to file
C.convertPyTree2File(a, 'out.cgns')