# CGNS Python tree basics

This tutorial explores the data handled by all Cassiopee functions,
called a *CGNS python tree* or shorter a *pytree*. This data is able to store the mesh coordinates, connectivities, fields and boundary conditions. It is conform to the CGNS standard document (https://cgns.org).

Get the sample cgns file: XXXX

## Reading a pytree from a file and basic tree node management

In this section, we will use the Converter module:

```python
import Converter.PyTree as C
import Converter.Internal as Internal
```

Here we read a sample file and display the tree to screen:

```python
# Read file
a = C.convertFile2PyTree("case.cgns")
# print pytree a to screen
Internal.printTree(a)
```

Here we get nodes of tree by their name or type:

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
is the node name, the second is the node value and is numpy, the third is all list
containing the child nodes and the last one is the node type. For example:

```python
import numpy
node = ['CoordinateX', numpy.array([0.,1.,2.], dtype=numpy.float64, order='F'), [], 'DataArray_t']
```

You can of course use numpy functions on numpy array (node[1]) as usual.

For example, change the X coordinate of the first point using numpy:
```python
# get numpy array of tree node
array = n1[1]
# change the X coordinate of first point
array[0] = 0.1
```

You can save the modified pytree to a file with:
```python
# save pytree to file
C.convertPyTree2File(a, 'out.cgns')
```

Internal offers other basic tree node management functions, such as node removal, node modification and node creation.

## Mesh topologies

The CGNS standard enables to store different mesh topotlogies: *structured* grids, *elements* grids and *polyedral* grids.

*Structured* grids are made of nixnjxnk points and dont need an explicit connectivity.
For example:

```python
import Generator.PyTree as G
# create a structured grid
a = G.cart((0,0,0), (1,1,1), (10,10,10))
Internal.printTree(a)
```
The internal numpy is then accessible by 3 indices.
```python
# get the CoordinateX node
n = Internal.getNodeFromName(a, 'CoordinateX')
print(n[1].shape)
```

Elements grids are of given element types (TRI, TETRA, HEXA, ...) and store a connectivity:
```python
# create an TETRA grid
a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
Internal.printTree(a)
```

This type of zone can also store different element types.
For example, create another grid with PENTA (prisms) and merge connectivities in a single multi-elements grid:

```python
# create a PENTA grid
b = G.cartPenta((0,0,9), (1,1,1), (10,10,5))
a = C.mergeConnectivity(a, b)
Internal.printTree(a)
```

Polyedral grids (also called NGON) enables faces of arbitrary number of points
and elements of arbitrary number of faces. Its is made of two connectivities,
the NGON described face node indices and the NFACE describes face indices.

```python
# create a NGON grid
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
Internal.printTree(a)
```

## Boundary conditions

Boundary conditions can be added topologically (from indices or range) or 
geometrically.

On structured grids, add a wall for "imin" window:

```python
# add a BC on a structured grid
a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._addBC2Zone(a, 'wall', 'BCWall', 'imin')
Internal.printTree(a)
```

On element or ngon grids, add a wall defined geometrically:

```python
a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
# define boundary geometry
b = G.cartHexa((0,0,0), (1,1,1), (10,10,1))
# add BC on unstructured grid
C._addBC2Zone(a, 'wall', 'BCWall', subzone=b)
Internal.printTree(a)
```

## Fields

In pytrees, you can add fields to any zones store in nodes or in centers. 
C.initVars enables to write formulas from other variables:

```python
# Init a node field
C._initVars(a, '{nodes:Density}=1.')
# Init a center field
C._initVars(a, '{centers:Pressure}=3*{centers:CoordinateX}')
Internal.printTree(a)
C.convertPyTree2File(a, 'out.cgns')
```

The "_" before initVars means that the function is applied "in place",
that is without duplication. Without "_", the function returns a reference copy:

```python
b = C.initVars(a, '{VelocityX}=0.')
```

You can also initialize a conservative field from a dimensional or non dimensional state:

```python
import Initiator.Adim as Adim
state = Adim.dim1(UInf=2.8, TInf=298., PInf=101325, LInf=12., alphaZ=1.)
C._initVars(a, 'centers:Density', state[0])
C._initVars(a, 'centers:MomentumX', state[1])
```