# - rotate (array) -
import Transform as T
import KCore.test as test

test.stdTestA(T.rotate, (0.,0.,0.), ((1.,0.,0.),(0,1,0),(0,0,1)),
              ((1,1,0), (1,-1,0), (0,0,1)))
