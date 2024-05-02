# - getFaceListOfEdge (PyTree) -
import Converter.Mpi as Cmpi
import OCC.PyTree as OCC
import KCore.test as test

edgeNo = 1
t = Cmpi.convertFile2PyTree("cube.cgns")
pos = OCC.getAllPos(t)

# Get face list from edgeNo
faceList = OCC.getFaceListOfEdge(t, pos, edgeNo)
test.testO(faceList,1)
