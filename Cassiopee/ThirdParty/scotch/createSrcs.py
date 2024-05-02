# Cree le fichier srcs contenant la liste des fichiers a compiler
import KCore.Dist as Dist

# cpp
cpp_srcs = Dist.getFilesOfExt('src/libscotch', ['.c'])

# f90
#f90_srcs = Dist.getFilesOfExt('Thermolib/BIBLIOTHEQUES', ['.f90'])

# tri pour trouver l'ordre de compilation (modules)
#f90_srcs = Dist.sortFileListByUse(f90_srcs)

f = open('srcs.py', 'w')
f.write('import KCore.Dist as Dist\n')
f.write('from KCore.config import *\n')

f.write('#==============================================================================\n')
f.write('# Fichiers c++\n')
f.write('#==============================================================================\n')
f.write('cpp_srcs = %s\n'%str(cpp_srcs))
f.write('#==============================================================================\n')

#f.write('# Fichiers fortran\n')
#f.write('#==============================================================================\n')
#f.write('f90_srcs = %s\n'%str(f90_srcs))
