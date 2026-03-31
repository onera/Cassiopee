# - adaptWindow2FacePointList (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

ni, nj, nk = 5, 5, 5

# inlet
win = [1, 1, 1, nj, 1, nk]
inletFPL = C.converter.window2FacePointList(*win, ni, nj, nk)
print(inletFPL)

# wall
win = [1, ni, 1, nj, 1, 1]
wallFPL = C.converter.window2FacePointList(*win, ni, nj, nk)
print(wallFPL)

# symm
win = [1, ni, nj, nj, 1, nk]
symmFPL = C.converter.window2FacePointList(*win, ni, nj, nk)
print(symmFPL)

a = G.cartHexa((0,0,0), (1,1,1), (ni, nj, nk))
C._addBC2Zone(a, 'inlet', 'BCInflow', faceList=inletFPL+1)
C._addBC2Zone(a, 'wall', 'BCWallViscous', faceList=wallFPL+1)
C._addBC2Zone(a, 'symm', 'BCSymmetry', faceList=symmFPL+1)
C.convertPyTree2File(a, 'out.cgns')
