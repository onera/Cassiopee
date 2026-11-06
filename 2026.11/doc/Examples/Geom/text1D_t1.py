# - text1D (array) -
import Geom as D
import KCore.test as test

# font text1
a1 = D.text1D("ABCDEFGHIJKLMNOPQRSTUVWXYZ", font='text1')
a2 = D.text1D("Aabcdefghijklmnopqrstuvwxyz", font='text1')
a3 = D.text1D("A0123456789+-=.,;:!()", font='text1')

a4 = D.text1D("ABCDEFGHIJKLMNOPQRSTUVWXYZ", font='text1', smooth=4, offset=1.)
a5 = D.text1D("A.1+2=3. 2-1=1.")
test.testA(a1, 1)

# font vera
a1 = D.text1D("ABCDEFGHIJKLMNOPQRSTUVWXYZ", font='vera')
a2 = D.text1D("Aabcdefghijklmnopqrstuvwxyz", font='vera')
a3 = D.text1D("A0123456789+-=.,;:!()", font='vera')
test.testA(a1, 2)

# Font chancery
a1 = D.text1D("ABCDEFGHIJKLMNOPQRSTUVWXYZ", font='chancery')
a2 = D.text1D("Aabcdefghijklmnopqrstuvwxyz", font='chancery')
a3 = D.text1D("A0123456789+-=.,;:!()", font='chancery')
test.testA(a1, 3)

# Font courier
a1 = D.text1D("ABCDEFGHIJKLMNOPQRSTUVWXYZ", font='courier')
a2 = D.text1D("Aabcdefghijklmnopqrstuvwxyz", font='courier')
a3 = D.text1D("A0123456789+-=.,;:!()", font='courier')
test.testA(a1, 4)

# Font nimbus
a1 = D.text1D("ABCDEFGHIJKLMNOPQRSTUVWXYZ", font='nimbus')
a2 = D.text1D("Aabcdefghijklmnopqrstuvwxyz", font='nimbus')
a3 = D.text1D("A0123456789+-=.,;:!()", font='nimbus')
test.testA(a1, 5)

test.writeCoverage(100)
