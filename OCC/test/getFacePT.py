# - getFace (PyTree) -
import Converter.Mpi as Cmpi
import OCC.PyTree as OCC

faceNo = 1
t = Cmpi.convertFile2PyTree("cube.cgns")
pos = OCC.getAllPos(t)
face1 = OCC.getFace(t, pos, faceNo)
