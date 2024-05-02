# - box (array) -
import Converter as C
import Transform as T
import Modeler.Models as Models

# Pure box
a = Models.box((0,3,0), (1,5,1))

# Box + chamfer
b = Models.box((0,0,0), (1,2,1), chamfer=0.1)

# Box2D
c = Models.box2D((2,2,0), (4,3,0), r=0.)
c = T.translate(c, (0,-2,0))

# Box2D + round
d = Models.box2D((2,2,0), (4,3,0), r=0.1)
d = T.translate(d, (0,0,0))

# Ellipse2D
e = Models.ellipse2D((2,2,0), (4,3,0))
e = T.translate(e, (0,2,0))

all = [a,b,c,d,e]
all = T.rotate(all, (0,0,0), (0,1,0), 90.)
C.convertArrays2File(all, 'out.plt')
