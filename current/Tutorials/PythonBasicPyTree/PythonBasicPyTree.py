
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

# create a node specifying name, data types and value
node = Internal.createNode('myNode', 'DataArray_t', value=[12.,14.,15.], children=[])
# create a node and attach it to another node as a child
child = Internal.createChild(node, 'myChild', 'DataArray_t', value=2.)
Internal.printTree(node)

# save pytree to file
C.convertPyTree2File(a, 'out.cgns')

# save pytree to tecplot format
C.convertPyTree2File(a, 'out.plt')
# save pytree to inria mesh format
C.convertPyTree2File(a, 'out.mesh', format='fmt_mesh')

import Generator.PyTree as G
# create a structured grid
a = G.cart((0,0,0), (1,1,1), (10,10,10))
Internal.printTree(a)
#> ['cart',array(shape=(3, 3),dtype='int32',order='F'),[2 sons],'Zone_t']
#>    |_['ZoneType',array('b'Structured'',dtype='|S1'),[0 son],'ZoneType_t']
#>    |_['GridCoordinates',None,[3 sons],'GridCoordinates_t']
#>        |_['CoordinateX',array(shape=(10, 10, 10),dtype='float64',order='F'),[0 son],'DataArray_t']
#>        |_['CoordinateY',array(shape=(10, 10, 10),dtype='float64',order='F'),[0 son],'DataArray_t']
#>        |_['CoordinateZ',array(shape=(10, 10, 10),dtype='float64',order='F'),[0 son],'DataArray_t']

# get the CoordinateX node
n = Internal.getNodeFromName(a, 'CoordinateX')
print(n[1].shape)
#> (10, 10, 10)

# create an TETRA zone
a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
Internal.printTree(a)
#> ['cartTetra',array(shape=(1, 3),dtype='int32',order='F'),[3 sons],'Zone_t']
#>    |_['ZoneType',array('b'Unstructured'',dtype='|S1'),[0 son],'ZoneType_t']
#>    |_['GridCoordinates',None,[3 sons],'GridCoordinates_t']
#>    |   |_['CoordinateX',array(shape=(1000,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>    |   |_['CoordinateY',array(shape=(1000,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>    |   |_['CoordinateZ',array(shape=(1000,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>    |_['GridElements',array(shape=(2,),dtype='int32',order='F'),[2 sons],'Elements_t']
#>        |_['ElementRange',array(shape=(2,),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>        |_['ElementConnectivity',array(shape=(14580,),dtype='int32',order='F'),[0 son],'DataArray_t']

# create a PENTA zone
b = G.cartPenta((0,0,9), (1,1,1), (10,10,5))
# merge zones a and b into a multi-element zone
a = C.mergeConnectivity(a, b)
Internal.printTree(a)
#> ['cartTetra',array(shape=(1, 3),dtype='int32',order='F'),[4 sons],'Zone_t']
#>    |_['ZoneType',array('b'Unstructured'',dtype='|S1'),[0 son],'ZoneType_t']
#>    |_['GridCoordinates',None,[3 sons],'GridCoordinates_t']
#>    |   |_['CoordinateX',array(shape=(1400,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>    |   |_['CoordinateY',array(shape=(1400,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>    |   |_['CoordinateZ',array(shape=(1400,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>    |_['GridElements',array(shape=(2,),dtype='int32',order='F'),[2 sons],'Elements_t']
#>    |   |_['ElementRange',array(shape=(2,),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>    |   |_['ElementConnectivity',array(shape=(14580,),dtype='int32',order='F'),[0 son],'DataArray_t']
#>    |_['GridElements-2',array(shape=(2,),dtype='int32',order='F'),[2 sons],'Elements_t']
#>        |_['ElementRange',array(shape=(2,),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>        |_['ElementConnectivity',array(shape=(3888,),dtype='int32',order='F'),[0 son],'DataArray_t']

# create a NGON grid
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
Internal.printTree(a)
#> ['cartNGon',array(shape=(1, 3),dtype='int32',order='F'),[4 sons],'Zone_t']
#>    |_['ZoneType',array('b'Unstructured'',dtype='|S1'),[0 son],'ZoneType_t']
#>    |_['GridCoordinates',None,[3 sons],'GridCoordinates_t']
#>    |   |_['CoordinateX',array(shape=(1000,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>    |   |_['CoordinateY',array(shape=(1000,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>    |   |_['CoordinateZ',array(shape=(1000,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>    |_['NGonElements',array(shape=(2,),dtype='int32',order='F'),[3 sons],'Elements_t']
#>    |   |_['ElementRange',array(shape=(2,),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>    |   |_['ElementConnectivity',array(shape=(12150,),dtype='int32',order='F'),[0 son],'DataArray_t']
#>    |   |_['FaceIndex',array(shape=(2430,),dtype='int32',order='F'),[0 son],'DataArray_t']
#>    |_['NFaceElements',array(shape=(2,),dtype='int32',order='F'),[3 sons],'Elements_t']
#>        |_['ElementRange',array(shape=(2,),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>        |_['ElementConnectivity',array(shape=(5103,),dtype='int32',order='F'),[0 son],'DataArray_t']
#>        |_['ElementIndex',array(shape=(729,),dtype='int32',order='F'),[0 son],'DataArray_t']

# add a BC on a structured grid
a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._addBC2Zone(a, 'wall', 'BCWall', 'imin')
Internal.printTree(a)
#> ['cart.0',array(shape=(3, 3),dtype='int32',order='F'),[3 sons],'Zone_t']
#>    |_['ZoneType',array('b'Structured'',dtype='|S1'),[0 son],'ZoneType_t']
#>    |_['GridCoordinates',None,[3 sons],'GridCoordinates_t']
#>    |   |_['CoordinateX',array(shape=(10, 10, 10),dtype='float64',order='F'),[0 son],'DataArray_t']
#>    |   |_['CoordinateY',array(shape=(10, 10, 10),dtype='float64',order='F'),[0 son],'DataArray_t']
#>    |   |_['CoordinateZ',array(shape=(10, 10, 10),dtype='float64',order='F'),[0 son],'DataArray_t']
#>    |_['ZoneBC',None,[1 son],'ZoneBC_t']
#>        |_['wall',array('b'BCWall'',dtype='|S1'),[1 son],'BC_t']
#>            |_['PointRange',array(shape=(3, 2),dtype='int32',order='F'),[0 son],'IndexRange_t']

a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
# define boundary geometry
b = G.cartHexa((0,0,0), (1,1,1), (10,10,1))
# add BC on unstructured grid
C._addBC2Zone(a, 'wall', 'BCWall', subzone=b)
Internal.printTree(a)
#> ['cartHexa.0',array(shape=(1, 3),dtype='int32',order='F'),[5 sons],'Zone_t']
#>    |_['ZoneType',array('b'Unstructured'',dtype='|S1'),[0 son],'ZoneType_t']
#>    |_['GridCoordinates',None,[3 sons],'GridCoordinates_t']
#>    |   |_['CoordinateX',array(shape=(1000,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>    |   |_['CoordinateY',array(shape=(1000,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>    |   |_['CoordinateZ',array(shape=(1000,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>    |_['GridElements',array(shape=(2,),dtype='int32',order='F'),[2 sons],'Elements_t']
#>    |   |_['ElementRange',array(shape=(2,),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>    |   |_['ElementConnectivity',array(shape=(5832,),dtype='int32',order='F'),[0 son],'DataArray_t']
#>    |_['cartHexa.1',array(shape=(2,),dtype='int32',order='F'),[2 sons],'Elements_t']
#>    |   |_['ElementRange',array(shape=(2,),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>    |   |_['ElementConnectivity',array(shape=(324,),dtype='int32',order='F'),[0 son],'DataArray_t']
#>    |_['ZoneBC',None,[1 son],'ZoneBC_t']
#>        |_['wall.0',array('b'BCWall'',dtype='|S1'),[1 son],'BC_t']
#>            |_['ElementRange',array(shape=(1, 2),dtype='int32',order='F'),[0 son],'IndexRange_t']

# Init a node field
C._initVars(a, '{nodes:Density}=1.')
# Init a center field
C._initVars(a, '{centers:Pressure}=3*{centers:CoordinateX}')
C.convertPyTree2File(a, 'out.cgns')

b = C.initVars(a, '{VelocityX}=0.')

import Initiator.Adim as Adim
state = Adim.dim1(UInf=2.8, TInf=298., PInf=101325, LInf=12., alphaZ=1.)
C._initVars(a, 'centers:Density', state[0])
C._initVars(a, 'centers:MomentumX', state[1])
