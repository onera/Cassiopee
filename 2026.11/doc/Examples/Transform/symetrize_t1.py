# - symmetrize (array) -
import Transform as T
import KCore.test as test

test.stdTestA(T.symmetrize, (0.,0.,0.), (1,0,0), (0,0,1))
