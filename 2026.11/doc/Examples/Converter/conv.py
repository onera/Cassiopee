# - convertFile2Arrays (binary tecplot) -
# - convertArrays2File (binary tecplot) -
# Read a file to a list of arrays
# Select some blocks in arrays
# Write the new arrays
import Converter as C

# Read file 'in.plt' to arrays
arrays = C.convertFile2Arrays("in.plt", "bin_tp")
n = len(arrays)
print('number of blocks in file: ', n)

# arrays is a list of n arrays.

# Write dimensions of first block : arrays[0]
print('dimensions: ', arrays[0][2], arrays[0][3], arrays[0][4])

# Write variable list of first block : arrays[0]
print('variables: ', arrays[0][0])

# Write coord. of first point of first block.
# ! : Note that index is not natural for Fortran programmers
print('x,y,z of first element: ',\
      arrays[0][1][0,0], arrays[0][1][1,0],  arrays[0][1][2,0])

# Select blocks 1 to 3 (first block is 0) and put them in a list
# named 'res'
res = arrays[1:4]

# Add block 5 at the end of list
res.append(arrays[5])

# Insert block 6 at beginning of list
res.insert(0, arrays[6])

# Save arrays to file 'out.plt' (bin_tp format)
C.convertArrays2File(res, 'out.plt')
