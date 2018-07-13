# - zipper (array) -
import Post as P
import KCore.test as test
import Converter as C

A = C.convertFile2Arrays("spoiler.plt")
a = P.zipper(A, ['overlapTol', 1.e-3])
test.testA([a], 1)
