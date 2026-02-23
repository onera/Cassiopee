# - text2D (array) -
import Geom as D
import KCore.test as test

myString = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz"
a = D.text2D(myString, smooth=0, offset=1., font='text1')
test.testA([a], 1)
a = D.text2D(myString, smooth=2, offset=1., font='text1')
test.testA([a], 2)
a = D.text2D(myString, smooth=0, offset=2., font='text1')
test.testA([a], 3)
a = D.text2D(myString, smooth=0, offset=1., font='vera')
test.testA([a], 4)

# Pas correct
a = D.text2D(myString, smooth=0, offset=2., font='chancery')
test.testA([a], 5)

a = D.text2D(myString, smooth=0, offset=1., font='courier')
test.testA([a], 6)
a = D.text2D(myString, smooth=0, offset=1., font='nimbus')
test.testA([a], 7)
test.writeCoverage(90)
