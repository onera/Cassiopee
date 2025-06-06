import pathlib
TIXI=False; TIGL=False

#==============================================================================
# Fichiers C++
#==============================================================================
cpp_srcs = []
if TIXI:
    srcs = [f.name for f in pathlib.Path("../../Modeler/tixi/src").glob("*.c")]
    for s in srcs: cpp_srcs.append("Modeler/tixi/src/"+s)
if TIGL:
    srcs = [f.name for f in pathlib.Path("../../Modeler/tigl/src").glob("*.cpp")]
    for s in srcs: cpp_srcs.append("Modeler/tigl/src/"+s)

    dirs = ['engine_nacelle', 'rotor', 'engine_pylon', 'structural_elements',
            'exports', 'system', 'api', 'fuelTanks', 'boolean_operations', 'fuselage',
            'common', 'generated', 'wing', 'configuration', 'geometry', 'contrib',
            'guide_curves', 'control_devices', 'imports', 'cpacs_other',
            'logging', 'ducts', 'math']
    for d in dirs:
        srcs = []
        for f in pathlib.Path("../../Modeler/tigl/src/"+d).glob("*.cpp"):
            srcs.append(f.name)
        for f in pathlib.Path("../../Modeler/tigl/src/"+d).glob("*.cxx"):
            srcs.append(f.name)
        for s in srcs: cpp_srcs.append("Modeler/tigl/src/%s/%s"%(d,s))
    print(cpp_srcs)

if TIGL:
    cpp_srcs += ['Modeler/CPACS/exportCAD.cpp']
else:
    cpp_srcs += ['Modeler/CPACS/exportCAD_stub.cpp']
