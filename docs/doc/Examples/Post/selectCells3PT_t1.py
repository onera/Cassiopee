# - selectCells3 (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Post.PyTree as P
import KCore.test as test

a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
C._initVars(a, '{centers:F}={centers:CoordinateX}')
C._initVars(a, '{centers:tag}=0.1*{centers:CoordinateX}+0.1*{centers:CoordinateY}')
b = P.selectCells3(a, 'centers:tag')
test.testT(b, 1)
