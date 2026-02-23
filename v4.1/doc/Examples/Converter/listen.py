# - listen (array) -
import Converter as C
import CPlot
sockets = C.createSockets()

while True:
    out = []
    for s in sockets:
        a = C.listen(s)
        if a is not None: out.append(a)
    if out != []: CPlot.display(out)
