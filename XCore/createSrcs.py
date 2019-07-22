# Cree le fichier srcs contenant la liste des fichiers a compiler
import KCore.Dist as Dist

# proper
cpp_srcs = []

# scotch
cpp_srcs1 = Dist.getFilesOfExt('XCore/scotch', ['.c'])

f = open('srcs.py', 'w')
f.write('import KCore.Dist as Dist\n')
f.write('from KCore.config import *\n')

f.write('#==============================================================================\n')
f.write('# Fichiers c++\n')
f.write('#==============================================================================\n')
f.write('cpp_srcs = %s\n'%str(cpp_srcs))
f.write('#==============================================================================\n')
f.write('cpp_srcs1 = %s\n'%str(cpp_srcs1))
f.write('#==============================================================================\n')
