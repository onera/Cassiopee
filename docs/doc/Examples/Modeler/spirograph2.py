# - spirograph (array) -
import Modeler.Models as Models
import CPlot

# k dans [0,1]: grandeur du cercle interieur
# l dans [0,1] : grandeur du pt A
# Ns: nbre de courbes
# N: nbre de pts sur chaque courbe
for i in range(100):
    for j in range(100):
        k = 0.01+i*0.01; l = j*0.01
        a = Models.spirograph(k, l, Ns=10, N=100)
        CPlot.display(a)
        CPlot.setState(message='k=%g, l=%g'%(k,l))
