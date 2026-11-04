#!/usr/bin/env python
import os
import sys
import KCore.test as test


# FILE_FORMATS: dictionary mapping each file extension to
#   - its format
#   - whether it is a PyTree (True) or an array (False).
FILE_FORMATS = {
    ".plt": ("bin_tp", False),
    ".dat": ("fmt_tp", False),
    ".tp": ("fmt_tp", False),
    ".v3d": ("bin_v3d", False),
    ".mesh": ("fmt_mesh", False),
    ".stl": ("bin_stl", False),
    ".fig": ("fmt_xfig", False),
    ".svg": ("fmt_svg", False),
    ".pov": ("fmt_pov", False),
    ".3ds": ("bin_3ds", False),
    ".obj": ("fmt_obj", False),
    ".gts": ("fmt_gts", False),
    ".cgns": ("bin_cgns", True),
    ".adf": ("bin_cgns", True),
    ".hdf": ("bin_hdf", True),
    ".pickle": ("bin_pickle", False),
}


def main():
    if len(sys.argv) < 3:
        print("Usage: kdiff file reffile [atol] [rtol]")
        sys.exit(1)
    # Parse arguments
    file1, file2 = sys.argv[1:3]
    if len(sys.argv) > 4:
        test.TOLERANCE = float(sys.argv[3])
    if len(sys.argv) == 5:
        test.RELTOLERANCE = float(sys.argv[4])

    for f in [file1, file2]:
        if not os.access(f, os.R_OK):
            print(f"Error: cannot read file '{f}'")
            sys.exit(1)
    ext = os.path.splitext(file1)[-1].lower()
    extInfo = FILE_FORMATS.get(ext, None)
    if extInfo is None:
        if ext[:4].lower() == ".ref":
            fmt, isPyTree = FILE_FORMATS[".pickle"]
        else:
            test.testO(file1, reference=file2)
            return
    else:
        fmt, isPyTree = extInfo

    if isPyTree:
        import Converter.PyTree as C
        t = C.convertFile2PyTree(file1, fmt)
        test.testT(t, reference=file2)
    else:
        import Converter
        a = Converter.convertFile2Arrays(file1, fmt)
        test.testA(a, reference=file2)
    return


if __name__ == "__main__":
    main()
