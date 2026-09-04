# - updateElementRange (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test
import numpy

# Create an HEXA with two BCs using subzones
a = G.cartHexa((0,0,0), (1,1,1), (5,5,5))
sz1 = G.cartHexa((0,0,0), (1,1,0), (5,5,1))
sz2 = G.cartHexa((0,0,4), (1,1,1), (5,5,1))
C._addBC2Zone(a, 'wall1', 'BCWall', subzone=sz1)
C._addBC2Zone(a, 'wall2', 'BCWall', subzone=sz2)
bc = Internal.getNodeFromName(a, 'wall1')
C._initBCDataSet(bc, 'Density', 0.)
bc = Internal.getNodeFromName(a, 'wall2')
C._initBCDataSet(bc, 'Density', 1.)
#C.convertPyTree2File(a, "out.cgns")

# Add a positive offset to the ElementRange of the second boundary Elements_t
# and second BC_t nodes such that sz1[1] != sz2[0] + 1
# Then call Internal._updateElementRange to fix it.
t = Internal.copyTree(a)
offset = 10
belt = Internal.getElementBoundaryNodes(t)[1]
ber = Internal.getNodeFromName1(belt, "ElementRange")[1]
ber += offset
bc = Internal.getNodeFromName(t, 'wall2')
ber = Internal.getNodeFromName1(bc, "ElementRange")[1]
ber += offset
Internal._updateElementRange(t)
test.testT(t, 1)
#C.convertPyTree2File(t, "out1.cgns")

# Delete boundary faces 0 and 5 from the last BC_t node and edit its
# corresponding boundary Elements_t node accordingly. Subsequently add a new
# boundary connectivity and BC. Indexing does not reauire a fix.
t = Internal.copyTree(a)
faceIndices2Delete = [0, 5]
belt = Internal.getElementBoundaryNodes(t)[1]
ber = Internal.getNodeFromName1(belt, "ElementRange")[1]
ber[1] -= len(faceIndices2Delete)
bec = Internal.getNodeFromName1(belt, "ElementConnectivity")
bec[1] = numpy.delete(bec[1].reshape(-1, 4), faceIndices2Delete, axis=0).ravel()

bc = Internal.getNodeFromName(t, 'wall2')
ber = Internal.getNodeFromName1(bc, "ElementRange")[1]
ber[0][1] -= len(faceIndices2Delete)
bcdata = Internal.getNodesFromType2(bc, "BCData_t")
for data in bcdata:
    datasets = Internal.getNodesFromType1(data, "DataArray_t")
    for ds in datasets:
        ds[1] = numpy.delete(ds[1], faceIndices2Delete)

sz3 = G.cartHexa((0,0,0), (0,1,1), (1,5,5))
C._addBC2Zone(t, 'wall3', 'BCWall', subzone=sz3)
test.testT(t, 2)
#C.convertPyTree2File(t, "out2.cgns")

# Delete boundary faces 0 and 5 from the BC_t node and edit its
# corresponding boundary Elements_t node accordingly. There is now a gap in
# ElementRange, call Internal._updateElementRange to fix it.
t = Internal.copyTree(a)
faceIndices2Delete = [0, 5]
belt = Internal.getElementBoundaryNodes(t)[0]
ber = Internal.getNodeFromName1(belt, "ElementRange")[1]
ber[1] -= len(faceIndices2Delete)
bec = Internal.getNodeFromName1(belt, "ElementConnectivity")
bec[1] = numpy.delete(bec[1].reshape(-1, 4), faceIndices2Delete, axis=0).ravel()

bc = Internal.getNodeFromName(t, 'wall1')
ber = Internal.getNodeFromName1(bc, "ElementRange")[1]
ber[0][1] -= len(faceIndices2Delete)
bcdata = Internal.getNodesFromType2(bc, "BCData_t")
for data in bcdata:
    datasets = Internal.getNodesFromType1(data, "DataArray_t")
    for ds in datasets:
        ds[1] = numpy.delete(ds[1], faceIndices2Delete)

Internal._updateElementRange(t)
test.testT(t, 3)
#C.convertPyTree2File(t, "out3.cgns")

# Add a negative offset to the ElementRange of the second boundary Elements_t
# and second BC_t nodes to create an overlap between ERs of sz1 and sz2 (such
# as one would get when adding faces to a boundary Elements_t).
# Then call Internal._updateElementRange to fix it.
t = Internal.copyTree(a)
offset = -10
belt = Internal.getElementBoundaryNodes(t)[1]
ber = Internal.getNodeFromName1(belt, "ElementRange")[1]
ber += offset
bc = Internal.getNodeFromName(t, 'wall2')
ber = Internal.getNodeFromName1(bc, "ElementRange")[1]
ber += offset
Internal._updateElementRange(t)
test.testT(t, 4)
#C.convertPyTree2File(t, "out4.cgns")
