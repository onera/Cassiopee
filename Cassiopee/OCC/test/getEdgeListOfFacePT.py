# - getEdgeListOfFace (PyTree) -
import Converter.Mpi as Cmpi
import OCC.PyTree as OCC

faceNo = 1
t = Cmpi.convertFile2PyTree("cube.cgns")
pos = OCC.getAllPos(t)
edgeList = OCC.getEdgeListOfFace(t, pos, faceNo)
