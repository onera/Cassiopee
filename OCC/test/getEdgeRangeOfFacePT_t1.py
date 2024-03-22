# - getEdgeRangeOfFace (PyTree) -
import Converter.Mpi as Cmpi
import OCC.PyTree as OCC
import KCore.test as test

faceNo = 1
t = Cmpi.convertFile2PyTree("cube.cgns")
pos = OCC.getAllPos(t)

# Return the ranges of edges of faceNo
r = OCC.getEdgeRangeOfFace(t, pos, faceNo)
test.testO(r,1)
