# - getEdgeRangeOfFace (PyTree) -
import Converter.Mpi as Cmpi
import OCC.PyTree as OCC

faceNo = 1
t = Cmpi.convertFile2PyTree("cube.cgns")
pos = OCC.getAllPos(t)
r = OCC.getEdgeRangeOfFace(t, pos, faceNo)
