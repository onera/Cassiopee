# Test KCore installation
# Remove . from PYTHONPATH
import sys, os
try:
    del sys.path[sys.path.index('')]
except: pass
try:
    del sys.path[sys.path.index(os.getcwd())]
except: pass
# Try to detect if colored output is supported
color = False
if sys.stdout.isatty(): color = True

# try import module
moduleName = 'KCore'
try:
    __import__(moduleName)
    if color:
        print("\033[32m%s correctly installed.\033[0m"%moduleName)
    else: print("%s correctly installed."%moduleName)
except Exception as inst:
    if color:
        print("\033[31mFAILED: %s\033[0m"%inst)
        print("\033[31mFAILED: %s badly installed.\033[0m"%moduleName)
    else:
        print("FAILED: %s"%inst)
        print("FAILED: %s badly installed."%moduleName)
