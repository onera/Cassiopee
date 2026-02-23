# - leave (array) -
import Converter as C
import Transform as T
import Modeler.Models as Models

# First leave
a = Models.leave(w=1., h=2.)

# Second leave
b = Models.leave(w=0.8, h=3.)
b = T.translate(b, (2,0,0))

C.convertArrays2File([a,b], 'out.plt')
