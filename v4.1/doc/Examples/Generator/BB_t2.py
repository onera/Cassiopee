# - BB (array) -
import Generator as G
import Transform as T
import Post as P
import KCore.test as test

# Rotated box
box = G.cartTetra((0,0,0), (1.,1.,1.), (2,2,2))
box = P.exteriorFaces(box)
box = T.rotate(box,(0,0,0),(1,0,0),+50.) # Rotation X
box = T.rotate(box,(0,0,0),(0,1,0),-20.) # Rotation Y
box = T.rotate(box,(0,0,0),(0,0,1),30.) # Rotation Z
OBB = G.BB(box,method='OBB',weighting=1)
test.testA([OBB],1)

# Aligned box
box = G.cartTetra((0,0,0), (1.,1.,1.), (2,2,2))
box = P.exteriorFaces(box)
OBB = G.BB(box,method='OBB',weighting=1)
test.testA([OBB],2)
