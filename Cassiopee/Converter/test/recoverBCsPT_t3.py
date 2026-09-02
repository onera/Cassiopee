# - recoverBCs (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test
import numpy

# recoverBCs with missing BCs
a = G.cartHexa((0,0,0),(1,1,1),(5,5,5))
sz1 = G.cartHexa((0,0,0),(1,1,0),(5,5,1))
sz2 = G.cartHexa((0,0,4),(1,1,1),(5,5,1))
C._addBC2Zone(a, 'wall1', 'BCWall', subzone=sz1)
C._addBC2Zone(a, 'wall2', 'BCWall', subzone=sz2)
#C.convertPyTree2File(a, "out.cgns")
t = Internal.copyTree(a)

# --- Geometric recoverBCs ---
# original unaltered zone
t = Internal.copyTree(a)
BCInfo = C.getBCs(t)
missingBCInfo = []
C._recoverBCsGeometric(t, BCInfo, removeBC=True,
                       missingBCInfo=missingBCInfo)
test.testO(missingBCInfo, 1)

# zone with modified boundary face connectivity and one corrupted BC
t = Internal.copyTree(a)
belt = Internal.getElementBoundaryNodes(t)[0]
bec = Internal.getNodeFromName(belt, "ElementConnectivity")[1]
bec[0] = 10
bec[-1] = 1

BCInfo = C.getBCs(t)
missingBCInfo = []
C._recoverBCsGeometric(t, BCInfo, removeBC=True,
                       missingBCInfo=missingBCInfo)
test.testO(missingBCInfo, 2)

# zone with modified boundary face connectivities and several corrupted BCs
t = Internal.copyTree(a)
belts = Internal.getElementBoundaryNodes(t)
for belt in belts:
    bec = Internal.getNodeFromName(belt, "ElementConnectivity")[1]
    bec[0] = 10
    bec[-1] = 1

BCInfo = C.getBCs(t)
missingBCInfo = []
C._recoverBCsGeometric(t, BCInfo, removeBC=True,
                       missingBCInfo=missingBCInfo)
test.testO(missingBCInfo, 3)

# --- Topologic recoverBCs ---
"""
# original unaltered zone
t = Internal.copyTree(a)
BCInfo = C.getBCs(t)
missingBCInfo = []
C._recoverBCsTopologic(t, BCInfo, removeBC=True,
                       missingBCInfo=missingBCInfo)
print(missingBCInfo)
test.testO(missingBCInfo, 1)

# zone with modified boundary face connectivity and one corrupted BC
t = Internal.copyTree(a)
belt = Internal.getElementBoundaryNodes(t)[0]
bec = Internal.getNodeFromName(belt, "ElementConnectivity")[1]
bec[0] = 10
bec[-1] = 1

BCInfo = C.getBCs(t)
missingBCInfo = []
C._recoverBCsTopologic(t, BCInfo, removeBC=True,
                       missingBCInfo=missingBCInfo)
print(missingBCInfo)
#test.testO(missingBCInfo, 2)

# zone with modified boundary face connectivities and several corrupted BCs
t = Internal.copyTree(a)
belts = Internal.getElementBoundaryNodes(t)
for belt in belts:
    bec = Internal.getNodeFromName(belt, "ElementConnectivity")[1]
    bec[0] = 10
    bec[-1] = 1

BCInfo = C.getBCs(t)
missingBCInfo = []
C._recoverBCsTopologic(t, BCInfo, removeBC=True,
                       missingBCInfo=missingBCInfo)
print(missingBCInfo)
#test.testO(missingBCInfo, 3)
"""
