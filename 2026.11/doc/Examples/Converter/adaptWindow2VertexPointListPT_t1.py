# - adaptWindow2VertexPointList (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

ni, nj, nk = 5, 5, 5

# inlet
win = [1, 1, 1, nj, 1, nk]
inletVPL = C.converter.window2VertexPointList(*win, ni, nj, nk)
test.testO(inletVPL, 1)

# wall
win = [1, ni, 1, nj, 1, 1]
wallVPL = C.converter.window2VertexPointList(*win, ni, nj, nk)
test.testO(wallVPL, 2)

# symm
win = [1, ni, nj, nj, 1, nk]
symmVPL = C.converter.window2VertexPointList(*win, ni, nj, nk)
test.testO(symmVPL, 3)

a = G.cartHexa((0,0,0), (1,1,1), (ni, nj, nk))
C._addBC2Zone(a, 'inlet', 'BCInflow', pointList=inletVPL+1)
C._addBC2Zone(a, 'wall', 'BCWallViscous', pointList=wallVPL+1)
C._addBC2Zone(a, 'symm', 'BCSymmetryPlane', pointList=symmVPL+1)
test.testT(a, 4)
