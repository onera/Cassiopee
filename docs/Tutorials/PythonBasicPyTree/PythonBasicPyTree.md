# CGNS Python tree basics

This tutorial explores the handling of data through Cassiopee functions,
called a *CGNS python tree* or shorter a *pytree*. The data can contain the mesh coordinates, connectivities, fields, and boundary conditions. It conforms to the latest CGNS standard (https://cgns.org).

## Reading a pytree from a file and basic tree node management

[Download the sample cgns file](case.cgns) 

In this section, we will use the Converter module:

```python
import Converter.PyTree as C
import Converter.Internal as Internal
```

Here we read a sample file and display the tree to the screen:

```python
# Read file
a = C.convertFile2PyTree("case.cgns")
# print pytree a to screen
Internal.printTree(a)
```

Here we get tree nodes from their name or type:

```python
# get a node by name (first met)
n1 = Internal.getNodeFromName(a, 'CoordinateX')
# get a node by type (first met)
n2 = Internal.getNodeFromType(a, 'Zone_t')
# get a list of nodes by type (all nodes)
n3 = Internal.getNodesFromType(a, 'FlowSolution_t')
# you can start search from any node
n4 = Internal.getNodeFromName(n3, 'Density')
```

Every tree node has the same structure and is a python list. The first element 
is the node name, the second is the node value and is in numpy, the third is all list
containing the child nodes, and the last one is the node type; e.g.:

```python
import numpy
node = ['CoordinateX', numpy.array([0.,1.,2.], dtype=numpy.float64, order='F'), [], 'DataArray_t']
```

You can of course use numpy functions on numpy arrays (node[1]) as one would normally do.

For example, change the X coordinate of the first point using numpy:
```python
# get numpy array of tree node
array = n1[1]
# change the X coordinate of first point
array[0] = 0.1
```

You can create nodes:
```python
# create a node specifying name, data types and value
node = Internal.createNode('myNode', 'DataArray_t', value=[12.,14.,15.], children=[])
# create a node and attach it to another node as a child
child = Internal.createChild(node, 'myChild', 'DataArray_t', value=2.)
Internal.printTree(node)
#> ['myNode',array(shape=(3,),dtype='float64',order='F'),[1 son],'DataArray_t']
#>    |_['myChild',array([2.0],dtype='float64'),[0 son],'DataArray_t']
```

Internal offers other basic tree node management functions, such as specific node creations, node removal, or node modifications.

You can save the modified pytree to a file with:
```python
# save pytree to file
C.convertPyTree2File(a, 'out.cgns')
```

You can also save it with different formats by changing the extension or by forcing the format explicitely. 
Neverthless, it should be noted that a lot of formats don't support all of the CGNS nodes, such as the boundary conditions. For example:

```python
# save pytree to tecplot format
C.convertPyTree2File(a, 'out.plt')
# save pytree to inria mesh format
C.convertPyTree2File(a, 'out.mesh', format='fmt_mesh')
```

## Mesh topologies

The CGNS standard enables the storage of different mesh topologies: *structured* grids, *elements* grids, and *polyedral* grids.

*Structured* grids are made of $n_i \times n_j \times n_k$ points and don't need an explicit connectivity.
Here, we use the Generator module to generate grids. For example:

```python
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
```
The internal numpy can then be accessible through the 3 indices (i,j,k).
```python
# get the CoordinateX node
n = Internal.getNodeFromName(a, 'CoordinateX')
print(n[1].shape)
#> (10, 10, 10)
```

Element grids are of given element types (TRI, TETRA, HEXA, ...) and store a connectivity:
```python
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
```

This type of zone can also store different element types.
For example, create another grid with PENTA (prisms) and merge the connectivities into a single multi-element grid:

```python
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
```

Polyhedral grids (also known as NGONs) enabls faces of an arbitrary number of points
and elements of an arbitrary number of faces. It contained two connectivities:
the NGON described face node indices and the NFACE describes face indices.

```python
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
```

## Boundary conditions

Boundary conditions can be added topologically (from indices or range) or 
geometrically.

On structured grids, add a wall at the "imin" window:

```python
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
```

On element or NGON grids, add a geometrically defined wall:

```python
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
```

## Fields

In pytrees, you can add fields to any zone and store them either at the cell vertices (known as nodes) or at the cell centers (known as centers). 
C.initVars support formula writing using other variables:

```python
# Init a node field
C._initVars(a, '{nodes:Density}=1.')
# Init a center field
C._initVars(a, '{centers:Pressure}=3*{centers:CoordinateX}')
C.convertPyTree2File(a, 'out.cgns')
```

The "\_" before initVars means that the function is applied "in place",
that is to say without duplication. Without "\_", the function returns a reference copy:

```python
b = C.initVars(a, '{VelocityX}=0.')
```

You can also initialize a conservative field from a dimensional or non-dimensional state:

```python
import Initiator.Adim as Adim
state = Adim.dim1(UInf=2.8, TInf=298., PInf=101325, LInf=12., alphaZ=1.)
C._initVars(a, 'centers:Density', state[0])
C._initVars(a, 'centers:MomentumX', state[1])
```
