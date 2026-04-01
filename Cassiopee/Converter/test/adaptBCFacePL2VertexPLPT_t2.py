# - adaptBCFacePL2VertexPL (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

# -- BE
a = G.cartHexa((0,0,0), (0.1,0.1,0.1), (5, 5, 5))
C._addBC2Zone(a, 'wall1', 'BCWallInviscid', faceList=[2 + 6*i for i in range(4)])
C._addBC2Zone(a, 'wall2', 'BCWallViscous', faceList=[74 + 6*i for i in range(4)])
C._initBCDataSet(a, '{var}=1.')

# Adapt face PL into vertex PL for all BC nodes, remove FaceCenter BCDatasets
b = Internal.copyTree(a)
b = Internal.adaptBCFacePL2VertexPL(b, remove=True)
test.testT(b, 1)

# Adapt face PL into vertex PL for all BC nodes with name wall* (in place), keep FaceCenter BCDatasets
b = Internal.copyTree(a)
bcs = Internal.getBCNodesFromName(b, 'wall*')
Internal._adaptBCFacePL2VertexPL(b, bcs=bcs, remove=False)
test.testT(b, 2)

# Adapt face PL into vertex PL for all BCWallViscous nodes (in place), keep FaceCenter BCDatasets
b = Internal.copyTree(a)
Internal._adaptBCFacePL2VertexPL(b, btype='BCWallViscous', remove=False)
test.testT(b, 3)

# -- 2D ME: tri-quad-tri
#                |
#               tri
"""a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (5,10,1))
b = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.2), (5,10,1))
c = G.cartTetra((0.8,0.,0.), (0.1,0.1,0.2), (5,10,1))
d = G.cartTetra((0.4,-0.9,0.), (0.1,0.1,0.2), (5,10,1))
a = C.mergeConnectivity([a, b, c, d], None)
t = C.newPyTree(['Base', a])
C._addBC2Zone(t, 'wall1', 'BCWallInviscid', faceList=[2 + 6*i for i in range(4)])
C._addBC2Zone(t, 'wall2', 'BCWallViscous', faceList=[74 + 6*i for i in range(4)])
C._initBCDataSet(t, '{var}=1.')
C.convertPyTree2File(t, 'out1.cgns')

# Adapt face PL into vertex PL for all BC nodes, remove FaceCenter BCDatasets
t = Internal.adaptBCFacePL2VertexPL(t, remove=True)
C.convertPyTree2File(t, 'out2.cgns')
#test.testT(t, 4)"""

# -- 3D ME: pyra - penta - hexa
"""a = G.cartPyra((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
b = G.cartPenta((0.4,0.,0.), (0.1,0.1,0.1), (5,5,5))
c = G.cartHexa((0.8,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.mergeConnectivity([a, b, c], None)
t = C.newPyTree(['Base', a])
C._addBC2Zone(t, 'wall1', 'BCWallInviscid', faceList=[2 + 6*i for i in range(4)])
#C._addBC2Zone(t, 'wall2', 'BCWallViscous', faceList=[74])
C._initBCDataSet(t, '{var}=1.')
C.convertPyTree2File(t, 'out1.cgns')

# Adapt face PL into vertex PL for all BC nodes, remove FaceCenter BCDatasets
t = Internal.adaptBCFacePL2VertexPL(t, remove=True)
C.convertPyTree2File(t, 'out2.cgns')
#test.testT(t, 5)"""
