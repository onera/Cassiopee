# - screw (array) -

import Converter as C
import Transform as T
import Modeler.Models as Models

# HEXA screw
a = Models.hexaScrew((0,0,0), r=1., h=0.4)

# HEXA screw with chamfer
b = Models.hexaScrew((2.5,0,0), r=1., h=0.4, chamfer=0.05)

# Rounded screw
c = Models.roundScrew((0,2.5,0), r=1., h=0.4)

# Rounded screw + drive
d = Models.roundScrew((2.5,2.5,0), r=1., h=0.4, drive=True)


all = [a,b,c,d]
all = T.rotate(all, (0,0,0), (0,1,0), 90.)
C.convertArrays2File(all, 'out.plt')