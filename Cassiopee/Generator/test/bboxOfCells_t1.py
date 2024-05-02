# - bboxOfCells (array) -
import Generator as G
import KCore.test as test

# test 2d structure
a1 = G.cart((0.,0.,0.),(0.1,0.1,1.),(20,20,1))
b = G.bboxOfCells(a1)
test.testA([b],1)

# test 3d structure
a2 = G.cart((0.,0.,0.),(0.1,0.1,1.),(20,20,20))
b = G.bboxOfCells(a2)
test.testA([b],2)

# test 2d tri
a3 = G.cartTetra((0.,0.,0.),(0.1,0.1,1.),(20,20,1))
b = G.bboxOfCells(a3)
test.testA([b],3)

# test 3d tetra
a4 = G.cartTetra((0.,0.,0.),(0.1,0.1,1.),(20,20,20))
b = G.bboxOfCells(a4)
test.testA([b],4)

# test 2d quad
a5 = G.cartHexa((0.,0.,0.),(0.1,0.1,1.),(20,20,1))
b = G.bboxOfCells(a5)
test.testA([b],5)

# test 3d hexa
a6 = G.cartHexa((0.,0.,0.),(0.1,0.1,1.),(20,20,20))
b = G.bboxOfCells(a6)
test.testA([b],6)

# test sur une liste
A = [a1,a2,a3]
B = G.bboxOfCells(A)
test.testA(B,7)
