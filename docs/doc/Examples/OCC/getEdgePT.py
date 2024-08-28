# - getEdge (PyTree) -
import Converter.Mpi as Cmpi
import OCC.PyTree as OCC

edgeNo = 1
t = Cmpi.convertFile2PyTree("cube.cgns")
pos = OCC.getAllPos(t)
edge1 = OCC.getEdge(t, pos, edgeNo)
