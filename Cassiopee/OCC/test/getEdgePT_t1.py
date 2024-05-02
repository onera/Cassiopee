# - getEdge (PyTree) -
import Converter.Mpi as Cmpi
import OCC.PyTree as OCC
import KCore.test as test

edgeNo = 1
t = Cmpi.convertFile2PyTree("cube.cgns")
pos = OCC.getAllPos(t)

# Retourne l'edge a partir de edgeNo (numero global CAD)
edge1 = OCC.getEdge(t, pos, edgeNo)
test.testT(edge1,1)
