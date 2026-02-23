# - getAllPos (PyTree) -
import Converter.Mpi as Cmpi
import OCC.PyTree as OCC

t = Cmpi.convertFile2PyTree("cube.cgns")
pos = OCC.getAllPos(t)
