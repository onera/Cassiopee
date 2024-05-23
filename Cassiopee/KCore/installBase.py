# This is the dictionary keeping track of installation.
# The key is the machine name. For each key a list is stored.
# [description,
# f77compiler, f90compiler, Cppcompiler, useOMP, static,
# additionalIncludePaths, additionalLibs, additionalLibPaths].
# Paths are list of strings. useOMP, static, useCuda are booleans.
# Others are strings.
installDict = {
###############################################################################
'DESKTOP...': [ 'Windows ubuntu',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/usr/include/mpi', '/usr/include/hdf5/serial'], # additionalIncludePaths
['gfortran', 'gomp', 'pthread'], # additionalLibs
['/usr/lib/x86_64-linux-gnu','/usr/lib/x86_64-linux-gnu/hdf5/serial/'], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'WDSNA81OZ': [ 'Machine de production win32 (Onera)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['c:/MinGW/include'], # additionalIncludePaths
['gfortran', 'gomp', 'pthread'], # additionalLibs
['c:/MinGW/lib', 'c:/Python27/libs', 'c:/MinGW/bin'], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'WDAAA728Z': [ 'Windows win64+msys2 (XJ-Onera)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
['-D__SHADERS__', '-isystem /d/juvigny/msys64/mingw64/include/python3.8', '-isystem /d/juvigny/msys64/mingw64/lib/python3.8/site-packages/numpy/core/include/numpy/'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['d:/juvigny/msys64/mingw64/include',"d:/juvigny/msys64/mingw64/include/OpenBLAS"], # additionalIncludePaths
['gfortran', 'gomp', 'pthread', 'openblas', 'psapi'], # additionalLibs
['d:/juvigny/msys64/mingw64/lib', 'd:/juvigny/msys64/mingw64/bin'], # additionalLibPaths
True, # useCuda
['-arch=sm_60'] # NvccAdditionalOptions
],
###############################################################################
'Aryen': [ 'Windows win64+msys2 (CB-Home)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
['-Wno-attributes', '-DSIMD=AVX2'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['c:/msys64/mingw64/include'], # additionalIncludePaths
['gomp', 'gfortran'], # additionalLibs
['c:/msys64/mingw64/lib'], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'msys64': [ 'Windows win64+msys2 (CB-Onera/Github)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
['-Wno-attributes', '-DSIMD=AVX2'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['d:/benoit/AppData/Local/msys2/mingw64/include', 'c:/Program Files (x86)/Microsoft SDKs/MPI/Include', 'd:/benoit/AppData/Local/msys2/mingw64/include/OpenBLAS'], # additionalIncludePaths
['gomp', 'gfortran'], # additionalLibs
['d:/benoit/AppData/Local/msys2/mingw64/lib', 'c:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64'], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'WDAAA878Z': [ 'Windows win64+msys2 (SL-Onera)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
['-Wno-attributes', '-fcommon'], # CppAdditionalOptions
[], # f77AdditionalOptions
False, # useOMP
False, # static
['c:/msys64/mingw64/include', 'c:/msys64/mingw64/include/OpenBLAS'], # additionalIncludePaths
['gomp', 'gfortran'], # additionalLibs
['c:/msys64/mingw64/lib'], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'node6.cluster': [ 'MacOSX (generic)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/usr/X11/include' ], # additionalIncludePaths
['python2.7', 'ifcore'], # additionalLibs
['/usr/X11/lib', '/System/Library/Frameworks/OpenGL.framework/Libraries/'], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'd1log1':[ 'Cluster HPC4B dev/val (Airbus)',
'ifort', # f77compiler
'ifort', # f90compiler
'icpc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/opt/hpmpi/include', '/opt/soft/cdtng/tools/portage/1.11/usr/include'], # additionalIncludePaths
['svml', 'irc', 'ifcore', 'iomp5'], # additionalLibs
['/opt/soft/cdtng/tools/intelcompiler/16.0/compiler/lib/intel64', '/opt/soft/cdtng/tools/portage/1.11/usr/lib', '/opt/hpmpi/lib/linux_amd64'], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'devt1n007.':[ 'Cluster HPC5 dev/val tou_b (Airbus)',
'ifort', # f77compiler
'ifort', # f90compiler
'icpc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/opt/mpi10/include', '/opt/soft/cdtng/tools/portage/1.11/usr/include'], # additionalIncludePaths
['ifcore', 'iomp5', 'svml', 'irc'], # additionalLibs
['/opt/soft/cdtng/tools/intelcompiler/2018/compilers_and_libraries_2018.5.274/linux/compiler/lib/intel64', '/opt/soft/cdtng/tools/portage/1.11/usr/lib', '/opt/mpi10/lib'], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'caefr0p...': [ 'Cluster GISEH (Airbus)',
'ifort', # f77compiler
'ifort', # f90compiler
'icpc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/opt/soft/cdtng/tools/portage/1.11/usr/include'], # additionalIncludePaths
['svml', 'irc', 'ifcore', 'iomp5'], # additionalLibs
['/opt/soft/cdtng/tools/portage/1.11/composerxe/lib/intel64', '/opt/soft/cdtng/tools/portage/1.11/usr/lib'],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'AA': [ 'Cluster AA',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=16'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/soft/intel-2017.5/compilers_and_libraries_2017.5.239/linux/mpi/include64','/home/mtaplaf/Environnement_Commun/Elsa/LIB_EXT/HDF5/hdf5-1.8.18-par/include'], # additionalIncludePaths
['mpi'], # additionalLibs
['/soft/intel-2017.5/compilers_and_libraries_2017.5.239/linux/mpi/lib64','/home/mtaplaf/Environnement_Commun/Elsa/LIB_EXT/HDF5/hdf5-1.8.18-par/lib'],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'AA2': [ 'Machine AA2',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=16'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/soft/intel-2017.5/compilers_and_libraries_2017.5.239/linux/mpi/include64','/home/mtaplaf/Environnement_Commun/Elsa/LIB_EXT/HDF5/hdf5-1.8.18-par/include'], # additionalIncludePaths
['mpi'], # additionalLibs
['/soft/intel-2017.5/compilers_and_libraries_2017.5.239/linux/mpi/lib64','/home/mtaplaf/Environnement_Commun/Elsa/LIB_EXT/HDF5/hdf5-1.8.18-par/lib'],  # additionalLibpaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'wfrontend1': [ 'Cluster Kairos (Safran)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
False, # useOMP
False, # static
['/appl/APPLI_SNECMA/HDF5/oper/1.8.11/include'], # additionalIncludePaths
['ifcore', 'svml', 'irc'], # additionalLibs
['/opt/intel/composer_xe_2013_sp1.0.080/lib/intel64', '/appl/APPLI_SNECMA/HDF5/oper/1.8.11/lib'],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'santafe': [ 'MacOSX - santafe (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/usr/X11/include' ], # additionalIncludePaths
['python2.7', 'ifcore'], # additionalLibs
['/usr/X11/lib', '/System/Library/Frameworks/OpenGL.framework/Libraries/'], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'daapuv': [ 'Machine DAAP (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/usr/local/hdf5-1.8.7/include'], # additionalIncludePaths
[], # additionalLibs
['/usr/local/hdf5-1.8.7/lib'],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'celeste': [ 'Grosse machine de post-traitement (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
['Xxf86vm'], # additionalLibs
[], # additionalLibPath
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'oneroa142': [ 'Machine dev (Onera)',
'/opt/intel/fc/9.1.036/bin/ifort', # f77compiler
'/opt/intel/fc/9.1.036/bin/ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'linux64': [ 'Production linux64 (generic)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/stck1/benoit/include'], # additionalIncludePaths
[], # additionalLibs
['/stck1/benoit/lib'],  # additionqlLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'eos814_r8': [ 'Poste grand calcul Onera-ld (Onera) avec Centos8',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
['-DCACHELINE=32','-DNB_SOCKET=1'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'eos8': [ 'Poste grand calcul Onera-ld (Onera) avec Centos8',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32','-DNB_SOCKET=1'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'eos...z': [ 'Poste grand calcul eosXXXz (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/usr/local/hdf5-gnu-1.8.8/include'], # additionalIncludePaths
[], # additionalLibs
['/usr/local/hdf5-gnu-1.8.8/lib'],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'eos...': [ 'Onera-eos (legacy-doit etre apres eosZ)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'ld...': [ 'Poste grand calcul Onera-ld (Onera)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
['-DCACHELINE=32','-DNB_SOCKET=1'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'tiamat': [ 'Machine de dev elsA (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=16','-DNB_SOCKET=1','-DCORE_PER_SOCK=6'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/home/benoit/x86_64t/include'], # additionalIncludePaths
[], # additionalLibs
['/home/benoit/x86_64t'],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'austri.onera': [ 'Cluster dev austri (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32','-DNB_SOCKET=2','-DCORE_PER_SOCK=12','-DSIMD=AVX2'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/home/benoit/aus/include'], # additionalIncludePaths
[], # additionalLibs
['/home/benoit/aus/lib'],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'rosetta-ws': [ 'Machine rosetta (Safran)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=16'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/softs/intel/compilers_and_libraries_2016.0.109/linux/mpi/intel64/include'], # additionalIncludePaths
[], # additionalLibs
['/softs/intel/compilers_and_libraries_2016.0.109/linux/mpi/intel64/lib'],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'rosetta-compute': [ 'Machine rosetta-calcul (Safran)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/softs/intel/compilers_and_libraries_2016.0.109/linux/mpi/intel64/include'], # additionalIncludePaths
[], # additionalLibs
['/softs/intel/compilers_and_libraries_2016.0.109/linux/mpi/intel64/lib'],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'westri': [ 'Machine westri-KNL (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=64','-DNB_SOCKET=1','-DCORE_PER_SOCK=64','-g','-DSIMD=MIC'], # CppAdditionalOptions
['-g'], # f77AdditionalOptions
True, # useOMP
False, # static
['/usr/local/hdf5-1.8.8-intel-16/include','/home/benoit/aus/include','/stck/nalferez/intel/parallel_studio_xe_2018/compilers_and_libraries_2018/linux/mpi/include64/','/stck/nalferez/intel/parallel_studio_xe_2018/vtune_amplifier_2018/include/','/stck/nalferez/intel/parallel_studio_xe_2018/advisor_2018/include/intel64'], # additionalIncludePaths
[], # additionalLibs
['/home/benoit/aus/lib','/stck/nalferez/intel/parallel_studio_xe_2018/compilers_and_libraries_2018/linux/mpi/lib64/','/stck/nalferez/intel/parallel_studio_xe_2018/compilers_and_libraries_2018/lib64/','/usr/local/hdf5-1.8.8-intel-16/lib/','/stck/nalferez/intel/parallel_studio_xe_2018/advisor_2018/lib64'],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'giulia': [ 'Machine dev elsA-ASO (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
["/tmp_opt/lib/hdf5-1.8.8-intel-16-impi/include",
"/usr/local/intel/studio/2016/compilers_and_libraries_2016.0.109/linux/mpi/include64"], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'mangrove': [ 'Machine avec acces GPU (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'moloch': [ 'Machine dev Cedre (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'cc-wdsna': [ 'Portable sous redhat (Onera)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'cephee': [ 'Cluster de dev Cassiopee (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/home/tools/local/x86_64a/include'], # additionalIncludePaths
[], # additionalLibs
['/home/tools/local/x86_64a/lib'],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'btmclx2': [ 'Cluster Turbomeca (Safran)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
['/usr/lib/gcc/x86_64-redhat-linux/4.1.2'],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'WDSNAXXX': [ '??',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'eos_pgi': [ 'Machine eos avec PGI',
'pgf90', # f77compiler
'pgf90', # f90compiler
'pgcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
["pgf90","pgf902","pgc","pgmath","pgf90_rpm1","rt","pgf90rtl","pgftnrtl"], # additionalLibs
["/d/juvigny/Logiciels/linux86-64/2018/lib"],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
##############################################################################
'visio': [ 'Machine de post gfx (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
['Xxf86vm'], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
##############################################################################
'visung': [ 'Machine de post gfx (Onera)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
##############################################################################
'visung_el8': [ 'Machine de post gfx (Onera)',
'ifx', # f77compiler
'ifx', # f90compiler
'icx', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
['Xxf86vm'], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'elmer': [ 'Machine de gros post gfx (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
['Xxf86vm'], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'WDSNA917Z': [ 'Machine de production win64 (Onera)',
'x86_64-w64-mingw32-gfortran', # f77compiler
'x86_64-w64-mingw32-gfortran', # f90compiler
'x86_64-w64-mingw32-gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
True, # static
['c:/TDM-GCC-64/include'], # additionalIncludePaths
['gfortran', 'gomp', 'quadmath'], # additionalLibs
['c:/TDM-GCC-64/lib', 'c:/Python2.7/libs'],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'fulvio': [ 'Machine post gfx legacy (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/usr/local/hdf5-intel-1.8.8/include'], # additionalIncludePaths
[], # additionalLibs
['/usr/local/hdf5-intel-1.8.8/lib'],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'cobalt': [ 'CCRT machine Cobalt',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'irene': [ 'TGCC machine Irene-skylake-partition',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'topaze': [ 'CCRT machine Topaze-milan-partition',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32','-DNB_SOCKET=2','-DCORE_PER_SOCK=64','-axCORE-AVX2','-mavx2'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionAllIbpaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'ouessant': [ 'Machine IDRIS IBM  POWER + NVIDIA P100)',
'pgf90', # f77compiler
'pgf90', # f90compiler
'pgcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
["pgf90","pgf902","pgc","pgmath","pgkomp",'omp',"pgf90_rpm1","rt","pgf90rtl"], # additionalLibs
["/opt/pgi/linuxpower/18.4/lib/"],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'jean-zay': [ 'Machine IDRIS intel + NVIDIA V100)',
'nvfortran', # f77compiler
'nvfortran', # f90compiler
'nvc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
True, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'ubuntu': [ 'Linux ubuntu (Onera)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'linux': [ 'Linux (generic)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'pdev': [ 'Machine Airbus (Airbus)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/opt/soft/cdtng/tools/portage/1.9/usr/include', '/opt/hpmpi/include'], # additionalIncludePaths
[], # additionalLibs
['/opt/soft/cdtng/tools/portage/1.9/usr/lib', '/opt/hpmpi/lib', '/opt/soft/cdtng/tools/intelcompiler/11.0/lib/intel64'],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'laura': [ 'Machine de dev acou (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/usr/local/hdf5/1.8.7/include'], # additionalIncludePaths
[], # additionalLibs
['/usr/local/lib64', '/usr/local/hdf5/1.8.7/lib','/tmp_opt/Python/2.7.3/icc-mpt/lib'],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'service': [ 'Cluster de calcul Stelvio (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
["-axAVX,SSE4.2"], # CppAdditionalOptions
["-axAVX,SSE4.2"], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
'r.i.n.': [ 'Cluster Stelvio-batch node (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
["-axAVX,SSE4.2"], # CppAdditionalOptions
["-axAVX,SSE4.2"], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'sator': [ 'Cluster de calcul Sator Broadwell (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32','-DNB_SOCKET=2','-DCORE_PER_SOCK=14','-Dvtune','-DSIMD=AVX2'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'sat_brw': [ 'Cluster de calcul Sator Broadwell (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32','-DNB_SOCKET=2','-DCORE_PER_SOCK=14','-Dvtune','-DSIMD=AVX2'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'sat_sky': [ 'Cluster de calcul Sator Skylake (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=64','-DNB_SOCKET=2','-DCORE_PER_SOCK=22','-Dvtune','-DSIMD=AVX2P512'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'sat_cas': [ 'Cluster de calcul Sator Cascadelake (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=64','-DNB_SOCKET=2','-DCORE_PER_SOCK=24','-Dvtune','-DSIMD=AVX2P512'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'spiro_sky': [ 'Machine dev Spiro (proc skylake)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=64','-DNB_SOCKET=2','-DCORE_PER_SOCK=22','-DSIMD=AVX512'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'spiro_amd': [ 'Machine dev Spiro (proc amd)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32','-DNB_SOCKET=2','-DCORE_PER_SOCK=32','-DSIMD=AVX2'],
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'spiro_arm': [ 'Machine dev Spiro (proc arm)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
['-DCACHELINE=32'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'spiro_acda': [ 'Machine dev Spiro (proc brwl)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'spiro_coda': [ 'Machine dev Spiro centos8 (Onera) (env. coda)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
['-DCACHELINE=32'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'spiro_el8': [ 'Machine dev Spiro centos8 (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32','-DNB_SOCKET=2','-DCORE_PER_SOCK=12','-DSIMD=AVX2'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'spiro_intel': [ 'Machine dev Spiro centos8 (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32','-DNB_SOCKET=2','-DCORE_PER_SOCK=12','-DSIMD=AVX2'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'spiro_gcc': [ 'Machine dev Spiro centos8 (Onera)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
['-DCACHELINE=32','-DNB_SOCKET=2','-DCORE_PER_SOCK=12','-DSIMD=AVX2'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'juno_gcc': [ 'Machine dev Juno rocky8 (Onera)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
['-DCACHELINE=32','-DNB_SOCKET=2','-DCORE_PER_SOCK=48','-DSIMD=AVX512'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'juno_coda': [ 'Machine dev Juno rocky8 (Onera) (env. coda)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
['-DCACHELINE=32'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'juno': [ 'Machine dev Juno rocky8 (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32','-DNB_SOCKET=2','-DCORE_PER_SOCK=48','-DSIMD=AVX512'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'pc_imad': [ 'pc imad local',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
#[],
['-DCACHELINE=32','-DNB_SOCKET=2','-DCORE_PER_SOCK=12','-DSIMD=AVX2'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'spiro_socle6': [ 'Machine dev Spiro centos8 - socle6 (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32','-DNB_SOCKET=2','-DCORE_PER_SOCK=12','-DSIMD=AVX2'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'spiro_pgi': [ 'Machine dev Spiro + compilos pgi (Onera)',
'nvfortran', # f77compiler
'nvfortran', # f90compiler
'nvc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
True, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'spiro': [ 'Machine dev Spiro (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32','-DNB_SOCKET=2','-DCORE_PER_SOCK=12','-DSIMD=AVX2'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'chi85bi': [ 'Cluster EDF (Onera)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'Raspail': [ 'Machine DTIM (Onera)',
'gfortran-7', # f77compiler
'gfortran-7', # f90compiler
'clang++-5.0', # Cppcompiler
#'g++-7', # Cppcompiler
['-pedantic', '-march=native', '-Wno-variadic-macros', '-Wno-long-long', '-g'], # CppAdditionalOptions
['-march=native', '-fdefault-real-8', '-fdefault-double-8'], # f77AdditionalOptions
True, # useOMP
False, # static
['/usr/include/hdf5/serial/'], # additionalIncludePaths
[], # additionalLibs
['/usr/lib/gcc/x86_64-linux-gnu/7',
 '/usr/lib/x86_64-linux-gnu/'],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'curie': [ 'Cluster Curie',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/usr/local/hdf5-1.8.8/include'], # additionalIncludePaths
[], # additionalLibs
['/usr/local/hdf5-1.8.8/lib'], # hdfPath
False, # useCuda
[] # NvccAdditionalOptions
],
 ##############################################################################
'madmax64': [ 'Cluster madmax DTIS (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
['/usr/local/intel/cluster_studio/2012_0_032/lib/intel64'], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'localhost.localdomain': [ 'Unknown',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'default': [ 'Default',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'stelvio_impi15': [ 'Cluster Stelvio Full intel (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
['/tmp_opt/lib/hdf5-1.8.8-intel-15-impi/include',
 # '/tmp_opt/lib/hdf5/1.8.17/15/impi/include',
 '/tmp_opt/intel/studio/2015/impi/5.0.3.048/intel64/include'], # additionalIncludePaths
['mpi'], # additionalLibs
['/tmp_opt/lib/hdf5-1.8.8-intel-15-impi/lib',
 '/tmp_opt/lib/hdf5/1.8.17/15/impi/lib',
 '/tmp_opt/intel/studio/2015/impi/5.0.3.048/intel64/lib'],  # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'dumbo': [ 'Grosse machine de post-traitement (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
['Xxf86vm'], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
],
###############################################################################
'xdaap': [ 'Xdaap (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
[], # additionalIncludePaths
[], # additionalLibs
[], # additionalLibPaths
False, # useCuda
[] # NvccAdditionalOptions
]
}
