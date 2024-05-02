# - reorder -
# Reordonne le bloc selectionne
# 's' sauvegarde le fichier
import Converter as C
import CPlot
import Transform as T
import time

a = C.convertFile2Arrays('fontaine.plt')
CPlot.display(a, mode=1, displayBB=0)

bool = 0
while (bool == 0):
    l = -1
    CPlot.display(a)
    while (l == -1):
        l = CPlot.getSelectedZone()
        s = CPlot.getKeyboard()
        if (s == "s"):
            C.convertArrays2File(a, 'out.plt')
            import sys; sys.exit();
        time.sleep(0.1)
    
    # Reorder suivant le type de zone
    z = a[l-1]
    if (len(z) == 5): # structure
        ni = z[2]; nj = z[3]; nk = z[4]
        if (ni == 1):
            a[l-1] = T.reorder(z, (1,2,-3))
        elif (nj == 1):
            a[l-1] = T.reorder(z, (1,2,-3))
        elif (nk == 1):
            a[l-1] = T.reorder(z, (-1,2,3))
        else:
            a[l-1] = T.reorder(z, (-1,2,3))
    else: # non structure
        a[l-1] = T.reorder(z, (-1,))

    CPlot.display(a)
