import pathlib
TIXI=False; TIGL=False

#==============================================================================
# Fichiers C++
#==============================================================================
cpp_srcs = []
if TIXI:
    cpp_srcs += [f.name for f in pathlib.Path("Tixi").glob("*.c")]
