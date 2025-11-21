# - addSquircle (array) -
import OCC

hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")
# add squircle
OCC.occ.addSuperEllipse(hook, (0,0,0), 1., 1., 4, 36, False)

# add super ellipse
OCC.occ.addSuperEllipse(hook, (0,3,0), 1., 1.5, 4, 36, False)

OCC.occ.writeCAD(hook, "out.step", "fmt_step")
