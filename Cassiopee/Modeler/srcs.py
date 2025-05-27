import pathlib
TIXI=False; TIGL=False

#==============================================================================
# Fichiers C++
#==============================================================================
cpp_srcs = []
if TIXI:
    srcs = [f.name for f in pathlib.Path("../../Modeler/tixi/src").glob("*.c")]
    for s in srcs: cpp_srcs.append("Modeler/tixi/src/"+s)
    print(cpp_srcs)