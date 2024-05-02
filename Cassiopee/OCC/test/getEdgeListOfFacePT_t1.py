# - getEdgeListOfFace (PyTree) -
import Converter.Mpi as Cmpi
import OCC.PyTree as OCC
import KCore.test as test

faceNo = 1
t = Cmpi.convertFile2PyTree("cube.cgns")
pos = OCC.getAllPos(t)

# Get edge list from faceNo
edgeList = OCC.getEdgeListOfFace(t, pos, faceNo)
test.testO(edgeList,1)
