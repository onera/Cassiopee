# - getEdgeVerticesOfFace (PyTree) -
import Converter.Mpi as Cmpi
import OCC.PyTree as OCC
import KCore.test as test

faceNo = 1
t = Cmpi.convertFile2PyTree("cube.cgns")
pos = OCC.getAllPos(t)

# Return the vertex indices for all edges of faceNo
inds = OCC.getEdgeVerticesOfFace(t, pos, faceNo)
test.testO(inds,1)
