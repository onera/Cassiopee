# - getFaceNoOppOfEdge (PyTree) -
import Converter.Mpi as Cmpi
import OCC.PyTree as OCC

faceNo = 1
edgeNo = 1
t = Cmpi.convertFile2PyTree("cube.cgns")
pos = OCC.getAllPos(t)
faceOppNo = OCC.getFaceNoOppOfEdge(t, pos, edgeNo, faceNo)
