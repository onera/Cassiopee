# Read a file 'output.plt' in centers
# Write each block to a separate file with bin_v3d format
import Converter as C

arrays = C.convertFile2Arrays("output.plt", "bin_tp")

i = 1
for ar in arrays:
    if i < 10:
        C.convertArrays2File([ar], "rep0"+repr(i)+".v3d", "bin_v3d")
    else:
        C.convertArrays2File([ar], "rep"+repr(i)+".v3d", "bin_v3d")
    i=i+1
