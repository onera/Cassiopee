# test driver
import OCC.Driver as D

# Create parameters implicitely
line1 = D.Line('line1', (0,0,0), (0,1,0))
line1.P[0].x.range = [-1,1]
#line1.print()
#line1.writeCAD("out.step")

line2 = D.Line('line2', (0,1,0), (1,1,0))
#line2.writeCAD("out.step")
#line2.P1 = line1.P2 # constraint

line3 = D.Line('line3', (1,1,0), (1,0,0))

line4 = D.Line('line4', (1,0,0), (0,0,0))

sketch = D.Sketch('sketch', [line1, line2, line3, line4])
sketch.writeCAD('out.step')
sketch.print()
