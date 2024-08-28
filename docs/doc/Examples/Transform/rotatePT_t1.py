# - rotate (pyTree) -
import Transform.PyTree as T
import KCore.test as test

test.stdTestT(T.rotate, (0.,0.,0.), (0.,0.,1.), 30.)
