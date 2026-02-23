# - getFields (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# STRUCT
z = G.cart((0.,0.,0.),(1.,1.,1.),(10,10,10))
res = C.getFields('GridCoordinates', z)
test.testA(res, 1)
z = C.initVars(z, 'F', 1)
z = C.initVars(z, 'centers:G', 2)
res = C.getFields('FlowSolution', z)
test.testA(res, 2)
res = C.getFields('FlowSolution#Centers', z)
test.testA(res, 3)

# BE
z = G.cartTetra((0.,0.,0.),(1.,1.,1.),(10,10,10))
z = C.initVars(z, 'F', 1)
z = C.initVars(z, 'centers:G', 2)
res = C.getFields('FlowSolution', z)
test.testA(res, 4)
res = C.getFields('FlowSolution#Centers', z)
test.testA(res, 5)

z = G.cartTetra((0.,0.,0.),(1.,1.,1.),(10,10,1))
z = C.initVars(z, 'F', 1)
z = C.initVars(z, 'centers:G', 2)
res = C.getFields('FlowSolution', z)
test.testA(res, 6)
res = C.getFields('FlowSolution#Centers', z)
test.testA(res, 7)

z = G.cartTetra((0.,0.,0.),(1.,1.,1.),(10,1,1))
z = C.initVars(z, 'F', 1)
z = C.initVars(z, 'centers:G', 2)
res = C.getFields('FlowSolution', z)
test.testA(res, 8)
res = C.getFields('FlowSolution#Centers', z)
test.testA(res, 9)

# NGON
z = G.cartNGon((0.,0.,0.),(1.,1.,1.),(10,10,10))
z = C.initVars(z, 'F', 1)
z = C.initVars(z, 'centers:G', 2)
res = C.getFields('FlowSolution', z)
test.testA(res, 10)
res = C.getFields('FlowSolution#Centers', z)
test.testA(res, 11)
