# - BB (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Post.PyTree as P
import KCore.test as test

# Rotated box
box = G.cartTetra((0,0,0), (1.,1.,1.), (2,2,2))
box = P.exteriorFaces(box)
box = T.rotate(box,(0,0,0),(1,0,0),+50.) # Rotation X
box = T.rotate(box,(0,0,0),(0,1,0),-20.) # Rotation Y
box = T.rotate(box,(0,0,0),(0,0,1),30.) # Rotation Z
OBB = G.BB(box,method='OBB',weighting=1)
test.testT(OBB,1)

# Aligned box
box = G.cartTetra((0,0,0), (1.,1.,1.), (2,2,2))
box = P.exteriorFaces(box)
OBB = G.BB(box,method='OBB',weighting=1)
test.testT(OBB,2)
