# - wall (array) -
import Converter as C
import Modeler.Models as Models
import Geom as D

a = Models.wall(D.line((0,0,0), (0,10,0)), 2., 1., 1, nlayers=10, chamfer=-1)

b = Models.wall(D.line((0,0,0), (10,20,0)), 2., 1., 1, nlayers=10, chamfer=0.1)

c = Models.wall(D.circle((0,-20,0), 10, tetas=0, tetae=120.), 2., 1., 1, nlayers=10,chamfer=0.1,shrink=0.9)

C.convertArrays2File([a,b,c], 'out.plt')
