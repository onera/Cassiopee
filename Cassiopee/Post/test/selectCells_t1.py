# - selectCells (array) -
import Post as P
import KCore.test as test

def F(x, y, z):
    if (x + y + z > 5.): return True
    else: return False
strict = 1
test.stdTestA(P.selectCells,F,[],['x','y','z'],strict)
