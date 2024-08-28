# - bboxOfCells (array) -
import Generator as G
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(20,20,20))
b = G.bboxOfCells(a)
print(b)
