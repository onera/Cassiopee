# generate a simple geometry
import Roms.Driver as D
import Roms.DB.DataBase as DataBase

T1 = D.Part("Part1")

# Create parameter
height = T1.Scalar("height")
height.range = [0.1, 2., 0.1]

# Create parametric geometry
P1 = T1.Point("P1", (0.,0.,0.))
P2 = T1.Point("P2", (1.,0.,0.))
P3 = T1.Point("P3", (1.,1.,0.))
T1.Eq(P3.y, height)
P4 = T1.Point("P4", (0.,1.,0.))
T1.Eq(P4.y, height)

line1 = T1.Line("line1", P1, P2)
line2 = T1.Line("line2", P2, P3)
line3 = T1.Line("line3", P3, P4)
line4 = T1.Line("line4", P4, P1)

sketch1 = T1.Sketch('sketch1', [line1, line2, line3, line4], h=[0.1,0.1,0.1])

# solve
T1.solve()

# create DB
db = DataBase.DataBase("BASE1", parameters=['height'])
T1.connect(db)
T1.instantiate({'height': 5.})
sketch1.writeCAD('out.step')
ref1 = sketch1.mesh()
db.registerReference(D.exportEdges(ref1), "ref1")

import CPlot, time
pt = T1.walkDOE()
while pt is not None:
    mesh = sketch1.rmesh(ref1)
    db.register("mesh", pt, "ref1", data=D.exportEdges(mesh))
    #CPlot.display(mesh)
    pt = T1.walkDOE()
    time.sleep(0.5)

# Fetch one
q = db.query('id = 1')
t = db.fetchTree(q)[0]
import Converter.PyTree as C
C.convertPyTree2File(t, 'out.cgns')

# Build full matrix from db
q = db.query()
A = db.fetchMatrix(q, variables=['dx','dy','dz'])
print(A.shape)