#!/usr/bin/env python

# You know how to write a script for converting files format:
# import Converter as C
# a = C.convertFile2Arrays("in.tp", "fmt_tp")
# C.convertArrays2File(a, "out.plt", "bin_tp")
#
# Now you can make it a shell command:

# First define a function equivalent to previous script:
def conv(fmttpFile, bintpFile):
    import Converter as C
    a = C.convertFile2Arrays(fmttpFile, "fmt_tp")
    C.convertArrays2File(a, bintpFile, "bin_tp")

# Then write a main function calling the previous function:
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("Two arguments are needed!")
    else:
        conv(sys.argv[1], sys.argv[2])
