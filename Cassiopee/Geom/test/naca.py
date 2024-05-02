# - naca (array) -
import Geom as D
import Converter as C

# Naca serie 4 defined by height
a = D.naca(12.)

# Naca serie 4 by name
b = D.naca('0012', N=301, sharpte=1)

# Naca serie 5 by name
c = D.naca('23012', N=301, sharpte=1)

# Naca serie 4 modified by name
d = D.naca('0008-45', N=301, sharpte=1)

C.convertArrays2File([a,b,c,d], 'out.plt')
