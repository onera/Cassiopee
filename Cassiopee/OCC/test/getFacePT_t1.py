# - getFace (PyTree) -
import Converter.Mpi as Cmpi
import OCC.PyTree as OCC
import KCore.test as test

faceNo = 1
t = Cmpi.convertFile2PyTree("cube.cgns")
pos = OCC.getAllPos(t)

# Retourne la face de faceNo (numero global CAD)
face1 = OCC.getFace(t, pos, faceNo)
test.testT(face1,1)
