# - square (array) -
from Generator.Shapes import square
import KCore.test as test

m = square((0,0), (1,1), 0.5, (10,10,5))
test.testA(m, 1)
