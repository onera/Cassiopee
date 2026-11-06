# - listen (pyTree) -
import Converter.PyTree as C
import CPlot.PyTree as CPlot

sockets = C.createSockets()
while True:
    out = []
    for s in sockets:
        a = C.listen(s)
        if a is not None: out.append(a)
    if out != []: CPlot.display(out)
