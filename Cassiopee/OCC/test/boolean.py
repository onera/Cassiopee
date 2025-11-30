# - boolean (array) -
import OCC

hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")

# first box
OCC.occ.addBox(hook, (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1))
# second box
OCC.occ.addBox(hook, (0.8,0.8,0.8), (1.2,0.8,0.8), (1.2,1.2,0.8), (0.8,1.2,0.8),
               (0.8,0.8,1.2), (1.2,0.8,1.2), (1.2,1.2,1.2), (0.8,1.2, 1.2))

# union
#OCC.occ.boolean(hook, [i for i in range(1,7)], [i for i in range(7,13)], 0, 1, 1)

# difference
#OCC.occ.boolean(hook, [i for i in range(1,7)], [i for i in range(7,13)], 1, 1, 1)

# intersection
OCC.occ.boolean(hook, [i for i in range(1,7)], [i for i in range(7,13)], 2, 1, 1)

OCC.occ.writeCAD(hook, "out.step", "fmt_step")
