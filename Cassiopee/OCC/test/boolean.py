# - boolean (array) -
import OCC

hook = OCC.createEmptyCAD()

# first box
OCC._addBox(hook, (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1))
# second box
OCC._addBox(hook, (0.8,0.8,0.8), (1.2,0.8,0.8), (1.2,1.2,0.8), (0.8,1.2,0.8),
            (0.8,0.8,1.2), (1.2,0.8,1.2), (1.2,1.2,1.2), (0.8,1.2, 1.2))

# union
#OCC._boolean(hook, [i for i in range(1,7)], [i for i in range(7,13)], 0, 1, 1)

# difference
#OCC._boolean(hook, [i for i in range(1,7)], [i for i in range(7,13)], 1, 1, 1)

# intersection
OCC._boolean(hook, [i for i in range(1,7)], [i for i in range(7,13)], 2, 1, 1)

OCC.writeCAD(hook, "out.step")
