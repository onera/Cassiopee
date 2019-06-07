# This is the dictionary keeping track of installation.
# The key is the machine name. For each key a list is stored.
# [description,
# f77compiler, f90compiler, Cppcompiler, useOMP, static, CPlotOffScreen,
# additionalIncludePaths, additionalLibs, additionalLibPaths].
# Paths are list of strings. useOMP, CPlotOffScreen are booleans.
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
False, # CPlotOffScreen
[], # additionalIncludePaths
['gfortran', 'gomp', 'pthread'], # additionalLibs
[] # additionalLibPaths
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
False, # CPlotOffScreen
['c:/MinGW/include'], # additionalIncludePaths
['gfortran', 'gomp', 'pthread'], # additionalLibs
['c:/MinGW/lib', 'c:/Python27/libs', 'c:/MinGW/bin'] # additionalLibPaths
],
###############################################################################
'WDAAA728Z': [ 'Windows win64+msys2 (XJ-Onera)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
['-isystem /d/juvigny/msys64/mingw64/include/python2.7', '-isystem /d/juvigny/msys64/mingw64/lib/python2.7/site-packages/numpy/core/include/numpy/'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
['d:/juvigny/msys64/mingw64/include',"c:/Program Files (x86)/Microsoft SDKs/MPI/Include", "/d/juvigny/msys64/mingw64/include/OpenBLAS"], # additionalIncludePaths
['gfortran', 'gomp', 'pthread', 'openblas', 'psapi'], # additionalLibs
['d:/juvigny/msys64/mingw64/lib', 'd:/juvigny/msys64/mingw64/bin',"c:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64"] # additionalLibPaths
],
###############################################################################
'WDAAA859Z': [ 'Windows win64+msys2 (CB-Onera)',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
['-Wno-attributes'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
['c:/msys64/mingw64/include', 'c:/Program Files (x86)/Microsoft SDKs/MPI/Include', 'c:/msys64/mingw64/include/OpenBLAS'], # additionalIncludePaths
['gomp', 'gfortran'], # additionalLibs
['c:/msys64/mingw64/lib', 'c:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64'] # additionalLibPaths
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
False, # CPlotOffScreen
['/usr/X11/include' ], # additionalIncludePaths
['python2.7', 'ifcore'], # additionalLibs
['/usr/X11/lib', '/System/Library/Frameworks/OpenGL.framework/Libraries/'], # additionalLibPaths
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
False, # CPlotOffScreen
['/opt/hpmpi/include', '/opt/soft/cdtng/tools/portage/1.11/usr/include'], # additionalIncludePaths
['svml', 'irc', 'ifcore', 'iomp5'], # additionalLibs
['/opt/soft/cdtng/tools/intelcompiler/16.0/compiler/lib/intel64', '/opt/soft/cdtng/tools/portage/1.11/usr/lib', '/opt/hpmpi/lib/linux_amd64'] # additionalLibPaths
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
False, # CPlotOffScreen
['/opt/mpi10/include', '/opt/soft/cdtng/tools/portage/1.11/usr/include'], # additionalIncludePaths
['svml', 'irc', 'ifcore', 'iomp5'], # additionalLibs
['/opt/soft/cdtng/tools/intelcompiler/2018/compilers_and_libraries_2018.5.274/linux/compiler/lib/intel64', '/opt/soft/cdtng/tools/portage/1.11/usr/lib', '/opt/mpi10/lib'] # additionalLibPaths
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
False, # CPlotOffScreen
['/opt/soft/cdtng/tools/portage/1.11/usr/include'], # additionalIncludePaths
['svml', 'irc', 'ifcore', 'iomp5'], # additionalLibs
['/opt/soft/cdtng/tools/portage/1.11/composerxe/lib/intel64', '/opt/soft/cdtng/tools/portage/1.11/usr/lib'] # additionalLibPaths
],
###############################################################################
'papin': [ 'Cluster MBDA (MBDA)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=16'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
True, # CPlotOffScreen
['/soft/intel-2017.5/compilers_and_libraries_2017.5.239/linux/mpi/include64','/home/mtaplaf/Environnement_Commun/Elsa/LIB_EXT/HDF5/hdf5-1.8.18-par/include'], # additionalIncludePaths
['mpi'], # additionalLibs
['/soft/intel-2017.5/compilers_and_libraries_2017.5.239/linux/mpi/lib64','/home/mtaplaf/Environnement_Commun/Elsa/LIB_EXT/HDF5/hdf5-1.8.18-par/lib'] # additionalLibPaths
],
###############################################################################
'eiffel': [ 'Machine MBDA2 (MBDA)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=16'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
True, # CPlotOffScreen
['/soft/intel-2017.5/compilers_and_libraries_2017.5.239/linux/mpi/include64','/home/mtaplaf/Environnement_Commun/Elsa/LIB_EXT/HDF5/hdf5-1.8.18-par/include'], # additionalIncludePaths
['mpi'], # additionalLibs
['/soft/intel-2017.5/compilers_and_libraries_2017.5.239/linux/mpi/lib64','/home/mtaplaf/Environnement_Commun/Elsa/LIB_EXT/HDF5/hdf5-1.8.18-par/lib'] # additionalLibPaths
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
False, # CPlotOffScreen
['/appl/APPLI_SNECMA/HDF5/oper/1.8.11/include'], # additionalIncludePaths
['ifcore', 'svml', 'irc'], # additionalLibs
['/opt/intel/composer_xe_2013_sp1.0.080/lib/intel64', '/appl/APPLI_SNECMA/HDF5/oper/1.8.11/lib'] # additionalLibPaths
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
False, # CPlotOffScreen
['/usr/X11/include' ], # additionalIncludePaths
['python2.7', 'ifcore'], # additionalLibs
['/usr/X11/lib', '/System/Library/Frameworks/OpenGL.framework/Libraries/'], # additionalLibPaths
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
True, # CPlotOffScreen
['/usr/local/hdf5-1.8.7/include'], # additionalIncludePaths
[], # additionalLibs
['/usr/local/hdf5-1.8.7/lib'] # additionalLibPaths
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
False, # CPlotOffScreen
[], # additionalIncludePaths
['Xxf86vm'], # additionalLibs
[] # additionalLibPath
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
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
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
False, # CPlotOffScreen
['/stck1/benoit/include'], # additionalIncludePaths
[], # additionalLibs
['/stck1/benoit/lib'] # additionalLibPaths
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
False, # CPlotOffScreen
['/usr/local/hdf5-gnu-1.8.8/include'], # additionalIncludePaths
[], # additionalLibs
['/usr/local/hdf5-gnu-1.8.8/lib'] # additionalLibPaths
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
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
],
###############################################################################
'ld...': [ 'Poste grand calcul Onera-ld (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32','-DNB_SOCKET=1'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
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
True, # CPlotOffScreen
['/home/benoit/x86_64t/include'], # additionalIncludePaths
[], # additionalLibs
['/home/benoit/x86_64t'] # additionalLibPaths
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
True, # CPlotOffScreen
['/home/benoit/aus/include'], # additionalIncludePaths
[], # additionalLibs
['/home/benoit/aus/lib'] # additionalLibPaths
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
True, # CPlotOffScreen
['/softs/intel/compilers_and_libraries_2016.0.109/linux/mpi/intel64/include'], # additionalIncludePaths
[], # additionalLibs
['/softs/intel/compilers_and_libraries_2016.0.109/linux/mpi/intel64/lib'] # additionalLibPaths
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
True, # CPlotOffScreen
['/softs/intel/compilers_and_libraries_2016.0.109/linux/mpi/intel64/include'], # additionalIncludePaths
[], # additionalLibs
['/softs/intel/compilers_and_libraries_2016.0.109/linux/mpi/intel64/lib'] # additionalLibPaths
],
###############################################################################
'westri': [ 'Machine westri-KNL (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=64','-DNB_SOCKET=1','-DCORE_PER_SOCK=64','-Dvtune','-g','-DSIMD=MIC'], # CppAdditionalOptions
['-g'], # f77AdditionalOptions
True, # useOMP
False, # static
True, # CPlotOffScreen
['/usr/local/hdf5-1.8.8-intel-16/include','/home/benoit/aus/include','/stck/nalferez/intel/parallel_studio_xe_2018/compilers_and_libraries_2018/linux/mpi/include64/','/stck/nalferez/intel/parallel_studio_xe_2018/vtune_amplifier_2018/include/','/stck/nalferez/intel/parallel_studio_xe_2018/advisor_2018/include/intel64'], # additionalIncludePaths
['ittnotify','advisor'], # additionalLibs
['/home/benoit/aus/lib','/stck/nalferez/intel/parallel_studio_xe_2018/compilers_and_libraries_2018/linux/mpi/lib64/','/stck/nalferez/intel/parallel_studio_xe_2018/compilers_and_libraries_2018/lib64/','/usr/local/hdf5-1.8.8-intel-16/lib/','/stck/nalferez/intel/parallel_studio_xe_2018/advisor_2018/lib64'] # additionalLibPaths
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
False, # CPlotOffScreen
["/tmp_opt/lib/hdf5-1.8.8-intel-16-impi/include",
"/usr/local/intel/studio/2016/compilers_and_libraries_2016.0.109/linux/mpi/include64"], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
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
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
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
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
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
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[], # additionalLibPaths
],
###############################################################################
'cephee': [ 'Cluser de dev Cassiopee (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
True, # CPlotOffScreen
['/home/tools/local/x86_64a/include'], # additionalIncludePaths
[], # additionalLibs
['/home/tools/local/x86_64a/lib'] # additionalLibPaths
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
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
['/usr/lib/gcc/x86_64-redhat-linux/4.1.2'] # additionalLibPaths
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
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
],
###############################################################################
'pgi': [ 'Machine eos avec PGI',
'pgf90', # f77compiler
'pgf90', # f90compiler
'pgcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
[], # additionalIncludePaths
["pgf90","pgf902","pgc","pgmath","pgf90_rpm1","rt","pgf90rtl","pgftnrtl"], # additionalLibs
["/d/juvigny/Logiciels/linux86-64/2018/lib"] # additionalLibPaths
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
False, # CPlotOffScreen
[], # additionalIncludePaths
['Xxf86vm'], # additionalLibs
[] # additionalLibPaths
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
False, # CPlotOffScreen
[], # additionalIncludePaths
['Xxf86vm'], # additionalLibs
[] # additionalLibPaths
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
False, # CPlotOffScreen
['c:/TDM-GCC-64/include'], # additionalIncludePaths
['gfortran', 'gomp', 'quadmath'], # additionalLibs
['c:/TDM-GCC-64/lib', 'c:/Python2.7/libs'] # additionalLibPaths
#['c:/TDM-GCC-64/lib', 'c:/Users/Adminstrateur/Anaconda2/libs'] # additionalLibPaths
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
False, # CPlotOffScreen
['/usr/local/hdf5-intel-1.8.8/include'], # additionalIncludePaths
[], # additionalLibs
['/usr/local/hdf5-intel-1.8.8/lib'] # additionalLibPaths
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
True, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
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
False, # CPlotOffScreen
[], # additionalIncludePaths
["pgf90","pgf902","pgc","pgmath","pgkomp",'omp',"pgf90_rpm1","rt","pgf90rtl"], # additionalLibs
["/opt/pgi/linuxpower/18.4/lib/"] # additionalLibPaths
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
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
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
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
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
False, # CPlotOffScreen
['/opt/soft/cdtng/tools/portage/1.9/usr/include', '/opt/hpmpi/include'], # additionalIncludePaths
[], # additionalLibs
['/opt/soft/cdtng/tools/portage/1.9/usr/lib', '/opt/hpmpi/lib', '/opt/soft/cdtng/tools/intelcompiler/11.0/lib/intel64'] # additionalLibPaths
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
False, # CPlotOffScreen
['/usr/local/hdf5/1.8.7/include'], # additionalIncludePaths
[], # additionalLibs
['/usr/local/lib64', '/usr/local/hdf5/1.8.7/lib','/tmp_opt/Python/2.7.3/icc-mpt/lib'] # additionalLibPaths
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
True, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
],
'r.i.n.': [ 'Cluster Stelvio-batch node (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
["-axAVX,SSE4.2"], # CppAdditionalOptions
["-axAVX,SSE4.2"], # f77AdditionalOptions
True, # useOMP
False, # static
True, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
],
###############################################################################
'sator': [ 'Cluster de calcul Sator (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32','-DNB_SOCKET=2','-DCORE_PER_SOCK=14','-Dvtune','-DSIMD=AVX2'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
True, # CPlotOffScreen
['/opt/tools/intel/studio/2017/compilers_and_libraries_2017.4.196/linux/mpi/intel64/include','/opt/tools/intel/studio/2017/vtune_amplifier_xe_2017.3.0.510739/include/','/opt/tools/intel/studio/2017/advisor_2017.1.3.510716/include/intel64'], # additionalIncludePaths
['ittnotify','advisor'], # additionalLibs
['/opt/tools/intel/studio/2017/vtune_amplifier_xe_2017.3.0.510739/lib64/','/opt/tools/intel/studio/2017/advisor_2017.1.3.510716/lib64','/opt/tools/lib/hdf5-1.8.17-intel-17/lib/'] # additionalLibPaths
#['/opt/tools/intel/studio/2017/vtune_amplifier_xe_2017.3.0.510739/lib64/','/opt/tools/intel/studio/2017/advisor_2017.1.3.510716/lib64'] # additionalLibPaths
],
###############################################################################
'sat_sky': [ 'Cluster de calcul Sator skylake(Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=64','-DNB_SOCKET=2','-DCORE_PER_SOCK=22','-Dvtune','-DSIMD=AVX512'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
True, # CPlotOffScreen
['/opt/tools/intel/studio/2017/compilers_and_libraries_2017.4.196/linux/mpi/intel64/include','/opt/tools/intel/studio/2017/vtune_amplifier_xe_2017.3.0.510739/include/','/opt/tools/intel/studio/2017/advisor_2017.1.3.510716/include/intel64'], # additionalIncludePaths
['ittnotify','advisor'], # additionalLibs
['/opt/tools/intel/studio/2017/vtune_amplifier_xe_2017.3.0.510739/lib64/','/opt/tools/intel/studio/2017/advisor_2017.1.3.510716/lib64','/opt/tools/lib/hdf5-1.8.17-intel-17/lib/'] # additionalLibPaths
#['/opt/tools/intel/studio/2017/vtune_amplifier_xe_2017.3.0.510739/lib64/','/opt/tools/intel/studio/2017/advisor_2017.1.3.510716/lib64'] # additionalLibPaths
],
###############################################################################
'spiro': [ 'Machine dev Spiro (Onera)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
True, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
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
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
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
False, # CPlotOffScreen
['/usr/include/hdf5/serial/'], # additionalIncludePaths
[], # additionalLibs
['/usr/lib/gcc/x86_64-linux-gnu/7',
 '/usr/lib/x86_64-linux-gnu/'] # additionalLibPaths
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
False, # CPlotOffScreen
['/usr/local/hdf5-1.8.8/include'], # additionalIncludePaths
[], # additionalLibs
['/usr/local/hdf5-1.8.8/lib'] # hdfPath
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
True, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
['/usr/local/intel/cluster_studio/2012_0_032/lib/intel64'], # additionalLibPaths
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
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
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
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
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
False, # CPlotOffScreen
['/tmp_opt/lib/hdf5-1.8.8-intel-15-impi/include',
 # '/tmp_opt/lib/hdf5/1.8.17/15/impi/include',
 '/tmp_opt/intel/studio/2015/impi/5.0.3.048/intel64/include'], # additionalIncludePaths
['mpi'], # additionalLibs
['/tmp_opt/lib/hdf5-1.8.8-intel-15-impi/lib', 
 '/tmp_opt/lib/hdf5/1.8.17/15/impi/lib', 
 '/tmp_opt/intel/studio/2015/impi/5.0.3.048/intel64/lib'] # additionalLibPaths
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
False, # CPlotOffScreen
[], # additionalIncludePaths
['Xxf86vm'], # additionalLibs
[] # additionalLibPath
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
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
]
}
