# - getAllPos (PyTree) -
import Converter.Mpi as Cmpi
import OCC.PyTree as OCC
import KCore.test as test

t = Cmpi.convertFile2PyTree("cube.cgns")

# Return the position of edges/faces in t
pos = OCC.getAllPos(t)
test.testO(pos,1)
