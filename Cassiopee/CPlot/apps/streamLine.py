# - streamline -
# trace des streamlines
import Converter as C
import CPlot
import Post as P
import time

a = C.convertFile2Arrays('outputn.plt')
b = C.convertFile2Arrays('output.plt')
CPlot.display(a, mode=0, displayBB=0, dim=2)

bool = 0
while (bool == 0):
    l = []
    CPlot.display(b)
    while (l == []):
        l = CPlot.getActivePoint()
        s = CPlot.getKeyboard()
        if (s == "s"):
            C.convertArrays2File(b, 'out.plt')
            import sys; sys.exit();
        time.sleep(0.1)
    stream = P.streamLine(a, (l[0],l[1],l[2])) 
    b += [stream]
    CPlot.display(b)
