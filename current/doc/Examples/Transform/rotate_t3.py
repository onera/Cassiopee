# - rotate (array) -
import Transform as T
import KCore.test as test

# test avec centre + 3 angles
test.stdTestA(T.rotate, (0.,0.,0.), (10.,20.,30.))
