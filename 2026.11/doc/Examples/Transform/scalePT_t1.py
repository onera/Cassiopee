# - scale (pyTree) -
import Transform.PyTree as T
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (11,11,11))

# scale in all directions
T._scale(a, factor=0.1)
test.testT(a, 1)

# scale in all directions with a reference point
T._scale(a, factor=0.1, X=(0,0,0))
test.testT(a, 2)

# scale with different factors following directions
T._scale(a, factor=(0.1,0.2,0.3))
test.testT(a, 3)
