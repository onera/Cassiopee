# - getNPtsOnEdgesOfFace (PyTree) -
import Converter.Mpi as Cmpi
import OCC.PyTree as OCC
import KCore.test as test

faceNo = 1
t = Cmpi.convertFile2PyTree("cube.cgns")
pos = OCC.getAllPos(t)

# Return the number of points on the edges of faceNo
nPts = OCC.getNPtsOnEdgesOfFace(t, pos, faceNo)
test.testO(nPts,1)
