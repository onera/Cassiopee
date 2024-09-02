# - box (array) -
import Transform as T
import Modeler.Models as Models
import KCore.test as test

# Pure box
a = Models.box((0,3,0), (1,5,1))

# Box + chamfer
b = Models.box((0,0,0), (1,2,1), chamfer=0.1)

# Box2D
c = Models.box2D((2,2,0), (5,4,0), r=0.)
c = T.rotate(c, (0,0,0), (0,1,0), 90.)

# Ellipse2D
d = Models.ellipse2D((2,2,0), (5,4,0))
d = T.rotate(d, (0,0,0), (0,1,0), 90.)

test.testA([a,b,c,d], 1)
