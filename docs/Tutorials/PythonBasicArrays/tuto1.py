import Converter as C
import Generator as G
import numpy

# Create a cartesian grid
a = G.cart( (0,0,0), (1,1,1), (10,11,12) )

# This is a structured Cassiopee array: ['x,y,z', n, ni, nj, nk]
# storing the variable string, the data as a numpy and the size of the
# grid
print(a)

# This is one way of checking type
if len(a) == 5: print('Structured')
else: print('Unstructured')

# Getting the numpy. This is a C indexed numpy array:
# n[0,ind] is x variable of index ind=i+j*ni+k*ni*nj
# n[1,ind] is y variable, and so on...
n = a[1]

# You can work on this numpy classicaly
n[2,:] = 2. # put all z to 2.

# Save the array in a file
# The bracket are here because this functions requires 'arrays' meaning
# a list of Cassiopee arrays.
C.convertArrays2File(a, 'out1.plt')

# Convert it as an hexa unstructured array
# Each function returns a copy of this original array
b = C.convertArray2Hexa(a)

# b is then a unstructured Cassiopee array: ['x,y,z', n, c, 'HEXA']
# storing the variable string, the data as a numpy and the connectivity
# as a numpy and the type of elements.
print(b)

C.convertArrays2File([a,b], 'out2.plt')
