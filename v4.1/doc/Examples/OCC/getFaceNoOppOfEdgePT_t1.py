# - getFaceNoOppOfEdge (PyTree) -
import Converter.PyTree as C
import OCC.PyTree as OCC
import KCore.test as test

faceNo = 1
edgeNo = 1
t = C.convertFile2PyTree("cube.cgns")
pos = OCC.getAllPos(t)

# Get the face opp of edgeNo belonging to faceNo
faceOppNo = OCC.getFaceNoOppOfEdge(t, pos, edgeNo, faceNo)
test.testO(faceOppNo,1)
