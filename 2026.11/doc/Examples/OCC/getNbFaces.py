# - getNbFaces (array) -
import OCC

hook = OCC.readCAD("cube.step", "fmt_step")
print(OCC.getNbFaces(hook))
#> 6
