# Cree le fichier srcs_paradigma contenant la liste des fichiers a compiler pour paradigma
import KCore.Dist as Dist

# paradigma
cpp_srcs = Dist.getFilesOfExt('XCore/paradigma', ['.c'])

f = open('srcs_paradigma.py', 'w')

f.write('#==============================================================================\n')
f.write('# Fichiers c++\n')
f.write('#==============================================================================\n')
f.write('cpp_srcs = %s\n'%str(cpp_srcs))
