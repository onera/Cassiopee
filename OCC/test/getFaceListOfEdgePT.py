# - getFaceListOfEdge (PyTree) -
import Converter.Mpi as Cmpi
import OCC.PyTree as OCC

edgeNo = 1
t = Cmpi.convertFile2PyTree("cube.cgns")
pos = OCC.getAllPos(t)
faceList = OCC.getFaceListOfEdge(t, pos, edgeNo)
