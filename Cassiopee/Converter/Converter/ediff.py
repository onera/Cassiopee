#!/usr/bin/env python
# Usage: ediff file1 file2
import sys
import Converter as C

def ediff__(file1, file2):
    a1 = C.convertFile2Arrays(file1)
    a2 = C.convertFile2Arrays(file2)
    ret = C.diffArrays(a1, a2)
    variables = C.getVarNames(ret)
    varsref = variables[0]
    for z in variables:
        for c, j in enumerate(z):
            if varsref[c] != j:
                print('Warning: ediff: variables are different in zones.')
    for i in varsref:
        L2 = C.normL2(ret, i)
        L0 = C.normL0(ret, i)
        if abs(L0) > 1.e-12:
            print('INFO : ', i, ', L2=', L2, ', L0=', L0)

def main():
    if len(sys.argv) < 3:
        print("Usage: ediff file1 file2")
        sys.exit(1)
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    ediff__(file1, file2)


if __name__ == "__main__":
    main()
