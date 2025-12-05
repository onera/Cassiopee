# test driver
import Roms.Driver as D
import Roms.DB.DataBase as DataBase

# Create parameter
height = D.Scalar("height", 1.)
height.range = [0.1, 2., 0.1]

P1 = D.Point("P1", (0.,0.,0.))
P2 = D.Point("P2", (1.,0.,0.))
P3 = D.Point("P3", (1.,height.v,0.))
D.Eq(P3.y.s, height.s)
P4 = D.Point("P4", (0.,height.v,0.))
D.Eq(P4.y.s, height.s)

line1 = D.Line("line1", P1, P2)
line2 = D.Line("line2", P2, P3)
line3 = D.Line("line3", P3, P4)
line4 = D.Line("line4", P4, P1)

sketch1 = D.Sketch('sketch1', [line1, line2, line3, line4])

# solve
D.DRIVER.solve2()

# create DB
db = DataBase.DataBase("BASE1", parameters=['height'])
D.DRIVER.connect(db)
D.DRIVER.instantiate({'height': 5.})
sketch1.writeCAD('out.step')
ref1 = sketch1.mesh(0.1, 0.1, 0.1)
db.registerReference(D.exportEdges(ref1), "ref1")

import CPlot, time
pt = D.DRIVER.walkDOE()
while pt is not None:
    mesh = sketch1.rmesh(ref1, 0.1, 0.1, 0.1)
    db.register("mesh", pt, "ref1", data=D.exportEdges(mesh))
    #CPlot.display(mesh)
    pt = D.DRIVER.walkDOE()
    time.sleep(0.5)

# Fetch one
q = db.query('id = 1')
t = db.fetchTree(q)[0]
import Converter.PyTree as C
C.convertPyTree2File(t, 'out.cgns')

# Create Model
q = db.query()
A = db.fetchMatrix(q, variables=['dx','dy','dz'])
print(A.shape)