# - spirograph (array) -
import Modeler.Models as Models
import Converter as C
import Transform as T

# k dans [0,1]: grandeur du cercle interieur
# l dans [0,1] : grandeur du pt A
# Ns: nbre de courbes
# N: nbre de pts sur chaque courbe
a = Models.spirograph(0.3, 0.6, Ns=10, N=100)
b = Models.spirograph(0.3, 0.99, Ns=10, N=100)
b = T.translate(b, (2,0,0))
c = Models.spirograph(0.18, 0.6, Ns=10, N=100)
c = T.translate(c, (0,2,0))
d = Models.spirograph(0.7, 0.2, Ns=10, N=100)
d = T.translate(d, (2,2,0))
C.convertArrays2File(a+b+c+d, 'out.plt')
