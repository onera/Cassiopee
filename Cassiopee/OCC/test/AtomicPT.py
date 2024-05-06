# - Atomic functions test -
import OCC
import Converter as C
import Generator as G

hook = OCC.occ.readCAD("cube.step", "fmt_step")

out = []
for i in range(6):
    # edges de la face i
    edges = OCC.occ.meshEdgesByFace(hook, i+1, 10, -1.)
    # edges dans espace uv
    edges = OCC.switch2UV(edges)
    # TFI dans espace uv
    a = G.TFI(edges)
    # evaluation sur la CAD
    o = OCC.occ.evalFace(hook, a, i+1)
    out.append(o)

C.convertArrays2File(out, 'mesh.plt')
