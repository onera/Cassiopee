
import Converter.PyTree as C
import Converter.Internal as Internal

# Read file
a = C.convertFile2PyTree("case.cgns")
# print pytree a to screen
Internal.printTree(a)

# get a node by name (first met)
n1 = Internal.getNodeFromName(a, 'CoordinateX')
# get a node by type (first met)
n2 = Internal.getNodeFromType(a, 'Zone_t')
# get a list of nodes by type (all nodes)
n3 = Internal.getNodesFromType(a, 'FlowSolution_t')
# you can start search from any node
n4 = Internal.getNodeFromName(n3, 'Density')

import numpy
node = ['CoordinateX', numpy.array([0.,1.,2.], dtype=numpy.float64, order='F'), [], 'DataArray_t']

# get numpy array of tree node
array = n1[1]
# change the X coordinate of first point
array[0] = 0.1

# save pytree to file
C.convertPyTree2File(a, 'out.cgns')

import Generator.PyTree as G
# create a structured grid
a = G.cart((0,0,0), (1,1,1), (10,10,10))
Internal.printTree(a)

# get the CoordinateX node
n = Internal.getNodeFromName(a, 'CoordinateX')
print(n[1].shape)

# create an TETRA grid
a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
Internal.printTree(a)

# create a PENTA grid
b = G.cartPenta((0,0,9), (1,1,1), (10,10,5))
a = C.mergeConnectivity(a, b)
Internal.printTree(a)

# create a NGON grid
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
Internal.printTree(a)

# add a BC on a structured grid
a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._addBC2Zone(a, 'wall', 'BCWall', 'imin')
Internal.printTree(a)

a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
# define boundary geometry
b = G.cartHexa((0,0,0), (1,1,1), (10,10,1))
# add BC on unstructured grid
C._addBC2Zone(a, 'wall', 'BCWall', subzone=b)
Internal.printTree(a)

# Init a node field
C._initVars(a, '{nodes:Density}=1.')
# Init a center field
C._initVars(a, '{centers:Pressure}=3*{centers:CoordinateX}')
Internal.printTree(a)
C.convertPyTree2File(a, 'out.cgns')

b = C.initVars(a, '{VelocityX}=0.')

import Initiator.Adim as Adim
state = Adim.dim1(UInf=2.8, TInf=298., PInf=101325, LInf=12., alphaZ=1.)
C._initVars(a, 'centers:Density', state[0])
C._initVars(a, 'centers:MomentumX', state[1])
