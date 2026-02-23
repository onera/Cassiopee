# - join (array) -
# test NACA0012
import Geom as D
import Transform as T
import Converter as C
import KCore.test as test

# Put a naca profile in a
a1 = D.naca(12., 5001)
a1 = C.initVars(a1, 'F', 2)

# Put a line in a2
a2 = D.line((1.,0.,0.), (20.,0.,0.), 5001)
a2 = C.initVars(a2, 'F', 3.)

# Join the two meshes
a1 = T.join(a1, a2)
# Add another line
a2 = D.line((1.,0.,0.), (20.,0.,0.), 5001)
a2 = C.initVars(a2, 'F', 3.)
a = T.join(a2, a1)
test.testA([a],1)
