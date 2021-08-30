import pathlib
import sys
this_source_dir = pathlib.Path(__file__).parent.resolve()

# do not print anything when importing
sys.stdout = open('/dev/null', 'w') # Something here that provides a write method.
import srcs
sys.stdout = sys.__stdout__

names = srcs.cpp_srcs + srcs.for_srcs + srcs.cpp_srcs2

for i in range(0,len(names)):
    names[i] = str(this_source_dir)+'/'+names[i];

print(';'.join(names))
#print("/scratchm/bberthou/travail/projects/KCore/KCore/Metis/balance.c")
