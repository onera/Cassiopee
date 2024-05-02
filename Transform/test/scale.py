# - scale (array) -
import Transform as T
import Generator as G
import Converter as C

a = G.cart((0,0,0), (1,1,1), (11,11,11))

# scale in all directions
a = T.scale(a, factor=0.1)

# scale in all directions with invariant point
a = T.scale(a, factor=0.1, X=(0,0,0))

# scale with different factors following directions
a = T.scale(a, factor=(0.1,0.2,0.3))

C.convertArrays2File(a, 'out.plt')
