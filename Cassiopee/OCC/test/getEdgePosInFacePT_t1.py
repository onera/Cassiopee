# - getEdgePosInFace (PyTree) -
import Converter.Mpi as Cmpi
import OCC.PyTree as OCC
import KCore.test as test

faceNo = 1
edgeNo = 1
t = Cmpi.convertFile2PyTree("cube.cgns")
pos = OCC.getAllPos(t)

# Return the position of edgeNo in faceNo
edgePos = OCC.getEdgePosInFace(t, pos, faceNo, edgeNo)
test.testO(edgePos,1)
