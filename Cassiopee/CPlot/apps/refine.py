# - refine -
# Raffine une cellule d'un maillage en clickant
# 's' sauvegarde le fichier
import Converter as C
import CPlot
import Post as P
import Generator as G
import time

a = [G.cart( (0,0,0), (1,1,1), (10,10,1) )]
a = C.convertArray2Tetra(a)
CPlot.display(a, mode=0, displayBB=0, dim=2)

bool = 0
while bool == 0:
    l = []
    CPlot.display(a)
    while l == []:
        l = CPlot.getActivePointIndex()
        nz = CPlot.getSelectedZone()
        s = CPlot.getKeyboard()
        if s == "s":
            C.convertArrays2File(a, 'out.plt')
            import sys; sys.exit();
        time.sleep(0.1)

    # Raffine
    z = a[nz]
    indic = C.array('indic', z[2].shape[1], 1, 1)
    indic = C.initVars(indic, 'indic', 0)
    C.setValue(indic, l[1], [1])
    a[nz] = P.refine(z, indic)

    CPlot.display(a)
