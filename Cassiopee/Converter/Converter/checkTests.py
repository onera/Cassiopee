# check unitary tests
# for every cassiopee functions
import sys, os, ast, subprocess

# Checks to perform
MISSINGTESTS = False
RUNTESTS = False
FUNCDOC = False
TESTHEADERS = False

# Modules to exclude
MODULEBLACKLIST = []
# python files to exclude
FILEBLACKLIST = ['__init__', 'cgnslib', 'cgnserrors', 'cgnskeywords', 'cgnstypes',
                 'cgnsutils', 'test',
                 'vera', 'text1', 'chancery', 'courier', 'nimbus']
# functions to exclude
FUNCTIONBLACKLIST = ['send']

# get the cassiopee source dir
def getCassiopeeSourceDir():
    try:
        # Check installPath
        import KCore.installPath
        cassiopeeIncDir = KCore.installPath.includePath
        cassiopeeIncDir = os.path.dirname(cassiopeeIncDir)
    except ImportError:
        raise SystemError("Error: KCore module is required to use this script.")
    return cassiopeeIncDir

# get the module list
def getModules():
    cassiopeeIncDir = getCassiopeeSourceDir()
    try: mods = os.listdir(cassiopeeIncDir)
    except: mods = []
    out = []
    for mod in mods:
        a = os.access(os.path.join(cassiopeeIncDir, mod, 'test'), os.F_OK)
        if a:
            if mod not in MODULEBLACKLIST:
                out.append(mod)
    return out

# get the module rst
def getModuleRST(module):
    cassiopeeIncDir = getCassiopeeSourceDir()
    docdir = os.path.join(cassiopeeIncDir, module, 'doc', 'source')
    docs = []
    if os.path.exists(docdir):
        files = os.listdir(docdir)
        for f in files:
            if f.endswith('.rst'): docs.append(os.path.join(docdir, f))
    doclines = []
    for doc in docs:
        with open(doc, "r", encoding="utf-8") as f:
            doclines += f.readlines()
    return doclines

# get the functions of module
# return a dict for each module python file containing
# (function name, function doc string)
def getFunctions(module):
    cassiopeeIncDir = getCassiopeeSourceDir()
    path = os.path.join(cassiopeeIncDir, module, module)
    files = os.listdir(path)
    modFuncs = {}
    for f in files: # python files in module

        if not f.endswith(".py"): continue
        name = f.replace(".py", "")
        # filter by file name
        if name in FILEBLACKLIST: continue

        modFuncs[name] = []
        file = os.path.join(path, f)
        with open(file, "r", encoding="utf-8") as f:
            tree = ast.parse(f.read())

        functions = []
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef):
                functions.append( (node.name, ast.get_docstring(node)) )
            if isinstance(node, ast.ClassDef):
                for sub in node.body:
                    if isinstance(sub, ast.FunctionDef):
                        functions.append( (f"{node.name}.{sub.name}", ast.get_docstring(sub)) )

        # filter by function names
        for fu in functions:
            if fu[0].endswith('__'): continue
            if fu[0][0] == '_': continue
            modFuncs[name].append(fu)

    return modFuncs

# IN: runTest: if true, run the tests of functionName
# IN: runHeader: if true, check test headers
# IN: runMissing: if true, check if tests are missing for functionName
# IN: runDoc: if true, check doc string and rst for functionName
def checkUnitaryTests(module, functionName, docString, docLines, runTest=False, runHeader=False,
                      runMissing=False, runDoc=False):

    if functionName in FUNCTIONBLACKLIST: return 0

    ret = 0
    cassiopeeIncDir = getCassiopeeSourceDir()
    path = os.path.join(cassiopeeIncDir, module, 'test')

    # check functionName docString
    if runDoc:
        if docString is None:
            print("Warning: %s docstring is missing."%(functionName))
            ret += 1
        if docLines is not None:
            found = False
            for l in docLines:
                if functionName in l:
                    found = True
                    break
            if not found:
                print("Warning: function %s is missing in rst."%(functionName))
                ret += 1

    # check name.py (doc test)
    name = functionName+'.py'
    file = os.path.join(path, name)
    a = os.path.exists(file)
    if not a:
        if runMissing:
            print("Warning: file %s is missing."%(name))
            ret += 1
    else:
        # check test
        if runTest:
            try:
                print("Running %s"%name)
                subprocess.check_output(["python", name], cwd=path, stderr=subprocess.STDOUT, text=True)
            except Exception:
                print("Error: test %s fails."%file)
                ret += 1

        # check header
        if runHeader:
            with open(file, "r", encoding="utf-8") as f:
                l = f.readline()
                ref = "# - %s (array) -"%functionName
                l = l[0:len(ref)]
                if l != ref:
                    print("Warning: bad header in file: %s."%file)
                    print("%s %s"%(ref, l))
                    ret += 1

    # check namePT.py (doc test)
    name = functionName+'PT.py'
    file = os.path.join(path, name)
    a = os.path.exists(file)
    if not a:
        if runMissing:
            print("Warning: file %s is missing."%(name))
            ret += 1
    else:
        # check test
        if runTest:
            try:
                subprocess.check_output(["python", name], cwd=path, stderr=subprocess.STDOUT, text=True)
            except Exception:
                print("Error: test %s fails."%file)
                ret += 1
        # check header
        if runHeader:
            with open(file, "r", encoding="utf-8") as f:
                l = f.readline()
                ref = "# - %s (pyTree) -"%functionName
                l = l[0:len(ref)]
                if l != ref:
                    print("Warning: bad header in file: %s."%file)
                    print("%s %s"%(ref, l))
                    ret += 1

    # check unitary test (array)
    name = functionName+'_t1.py'
    file = os.path.join(path, name)
    a = os.path.exists(file)
    if not a:
        if runMissing:
            print("Warning: file %s is missing."%(name))
            ret += 1
    else:
        # check header
        if runHeader:
            with open(file, "r", encoding="utf-8") as f:
                l = f.readline()
                ref = "# - %s (array) -"%functionName
                l = l[0:len(ref)]
                if l != ref:
                    print("Warning: bad header in file: %s."%file)
                    print("%s %s"%(ref, l))
                    ret += 1

    # check unitary test (pyTree)
    name = functionName+'PT_t1.py'
    file = os.path.join(path, name)
    a = os.path.exists(file)
    if not a:
        if runMissing:
            print("Warning: file %s is missing."%(name))
            ret += 1
    else:
        # check header
        if runHeader:
            with open(file, "r", encoding="utf-8") as f:
                l = f.readline()
                ref = "# - %s (pyTree) -"%functionName
                l = l[0:len(ref)]
                if l != ref:
                    print("Warning: bad header in file: %s."%file)
                    print("%s %s"%(ref, l))
                    ret += 1
    return ret

#==============================================================================
if __name__ == "__main__":
    argv = sys.argv
    if len(argv) == 1: # pas de module specifie
        modules = getModules()
    else:
        modules = argv[1:]
        cassiopeeIncDir = getCassiopeeSourceDir()
        for module in modules:
            path = os.path.join(cassiopeeIncDir, module, 'test')
            a = os.path.exists(path)
            if not a:
                raise ValueError("Module %s not found."%module)

    errorTot = 0
    for module in modules: # for each module
        errors = 0
        docLines = getModuleRST(module)
        funcs = getFunctions(module)
        fileNames = funcs.keys()
        for file in fileNames: # for each python file
            functionNames = funcs[file]
            for functionName in functionNames:
                e = checkUnitaryTests(module, functionName[0], functionName[1], docLines,
                                      runMissing=MISSINGTESTS, runHeader=TESTHEADERS,
                                      runTest=RUNTESTS, runDoc=FUNCDOC)
                errors += e
        if errors > 0:
            print("%s has %d python files."%(module, len(fileNames)))
            for file in fileNames:
                print("     %s/%s has %d functions."%(module, file, len(functionNames)))
            print("%s has %d errors."%(module, errors))
        errorTot += errors
    print("Total of %s errors."%errorTot)
