# - deformPoint (array) -
import Transform as T
import KCore.test as test

test.stdTestA(T.deformPoint, (0.,0.,0.), (0.1,0.1,0.), 0.5, 0.4)
