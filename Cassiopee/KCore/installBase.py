# This is the dictionary keeping track of installation.
# The key is the machine name. For each key a list is stored.
# [description,
# f77compiler, f90compiler, Cppcompiler, useOMP, static,
# additionalIncludePaths, additionalLibs, additionalLibPaths].
# Paths are list of strings. useOMP, static, useCuda are booleans.
# Others are strings.
try:
    from installBaseUser import installDict as installDictUser
except ImportError:
    try:
        from . import installBaseUser
        installDictUser = installBaseUser.installDict
    except:
        installDictUser = {}

installDict = {
    **installDictUser,
    ###############################################################################
    'WDAAA728Z': [
        'Windows win64+msys2 (XJ-ONERA)',
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
    'Aryen': [
        'Windows win64+msys2 (CB-Home)',
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
    'win64': [
        'Windows win64+msys2',
        'gfortran', # f77compiler
        'gfortran', # f90compiler
        'gcc', # Cppcompiler
        ['-Wno-attributes'], # CppAdditionalOptions
        [], # f77AdditionalOptions
        True, # useOMP
        False, # static
        [], # additionalIncludePaths
        ['gomp', 'gfortran'], # additionalLibs
        [], # additionalLibPaths
        False, # useCuda
        [] # NvccAdditionalOptions
    ],
    ###############################################################################
    'WDAAA161Z': [
        'Windows win64+msys2 (CB-ONERA)',
        'gfortran', # f77compiler
        'gfortran', # f90compiler
        'gcc', # Cppcompiler
        ['-Wno-attributes', '-fcommon'], # CppAdditionalOptions
        [], # f77AdditionalOptions
        True, # useOMP
        False, # static
        ['c:/Users/benoit/msys64/mingw64/include', 'c:/Users/benoit/msys64/mingw64/include/openblas'], # additionalIncludePaths
        ['gomp', 'gfortran'], # additionalLibs
        ['c:/Users/benoit/msys64/mingw64/lib'], # additionalLibPaths
        False, # useCuda
        [] # NvccAdditionalOptions
    ],
    ###############################################################################
    'node6.cluster': [
        'MacOSX (generic)',
        'ifort', # f77compiler
        'ifort', # f90compiler
        'icc', # Cppcompiler
        [], # CppAdditionalOptions
        [], # f77AdditionalOptions
        True, # useOMP
        False, # static
        ['/usr/X11/include'], # additionalIncludePaths
        ['python2.7', 'ifcore'], # additionalLibs
        ['/usr/X11/lib', '/System/Library/Frameworks/OpenGL.framework/Libraries/'], # additionalLibPaths
        False, # useCuda
        [] # NvccAdditionalOptions
    ],
    ###############################################################################
    'd1log1': [
        'Cluster HPC4B dev/val (Airbus)',
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
    'devt1n007.': [
        'Cluster HPC5 dev/val tou_b (Airbus)',
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
    'AA': [
        'Cluster AA',
        'ifort', # f77compiler
        'ifort', # f90compiler
        'icc', # Cppcompiler
        ['-DCACHELINE=16'], # CppAdditionalOptions
        [], # f77AdditionalOptions
        True, # useOMP
        False, # static
        ['/soft/intel-2017.5/compilers_and_libraries_2017.5.239/linux/mpi/include64', '/home/mtaplaf/Environnement_Commun/Elsa/LIB_EXT/HDF5/hdf5-1.8.18-par/include'], # additionalIncludePaths
        ['mpi'], # additionalLibs
        ['/soft/intel-2017.5/compilers_and_libraries_2017.5.239/linux/mpi/lib64', '/home/mtaplaf/Environnement_Commun/Elsa/LIB_EXT/HDF5/hdf5-1.8.18-par/lib'], # additionalLibPaths
        False, # useCuda
        [] # NvccAdditionalOptions
    ],
    ###############################################################################
    'AA2': [
        'Machine AA2',
        'ifort', # f77compiler
        'ifort', # f90compiler
        'icc', # Cppcompiler
        ['-DCACHELINE=16'], # CppAdditionalOptions
        [], # f77AdditionalOptions
        True, # useOMP
        False, # static
        ['/soft/intel-2017.5/compilers_and_libraries_2017.5.239/linux/mpi/include64', '/home/mtaplaf/Environnement_Commun/Elsa/LIB_EXT/HDF5/hdf5-1.8.18-par/include'], # additionalIncludePaths
        ['mpi'], # additionalLibs
        ['/soft/intel-2017.5/compilers_and_libraries_2017.5.239/linux/mpi/lib64', '/home/mtaplaf/Environnement_Commun/Elsa/LIB_EXT/HDF5/hdf5-1.8.18-par/lib'], # additionalLibPaths
        False, # useCuda
        [] # NvccAdditionalOptions
    ],
    ###############################################################################
    'wfrontend1': [
        'Cluster Kairos (Safran)',
        'ifort', # f77compiler
        'ifort', # f90compiler
        'icc', # Cppcompiler
        [], # CppAdditionalOptions
        [], # f77AdditionalOptions
        False, # useOMP
        False, # static
        ['/appl/APPLI_SNECMA/HDF5/oper/1.8.11/include'], # additionalIncludePaths
        ['ifcore', 'svml', 'irc'], # additionalLibs
        ['/opt/intel/composer_xe_2013_sp1.0.080/lib/intel64', '/appl/APPLI_SNECMA/HDF5/oper/1.8.11/lib'], # additionalLibPaths
        False, # useCuda
        [] # NvccAdditionalOptions
    ],
    ###############################################################################
    'celeste': [
        'Grosse machine de post-traitement (ONERA)',
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
    'linux64': [
        'Production linux64 (generic)',
        'gfortran', # f77compiler
        'gfortran', # f90compiler
        'gcc', # Cppcompiler
        [], # CppAdditionalOptions
        [], # f77AdditionalOptions
        True, # useOMP
        False, # static
        ['/stck1/benoit/include'], # additionalIncludePaths
        [], # additionalLibs
        ['/stck1/benoit/lib'], # additionalLibPaths
        False, # useCuda
        [] # NvccAdditionalOptions
    ],
    ###############################################################################
    'ld_spack1': [
        'Poste grand calcul Onera-ld (ONERA)',
        'gfortran', # f77compiler
        'gfortran', # f90compiler
        'gcc', # Cppcompiler
        ['-DCACHELINE=32', '-DNB_SOCKET=1'], # CppAdditionalOptions
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
    'ld_spack2': [
        'Poste grand calcul Onera-ld (ONERA)',
        'ifort', # f77compiler
        'ifort', # f90compiler
        'icc', # Cppcompiler
        ['-DCACHELINE=32', '-DNB_SOCKET=1'], # CppAdditionalOptions
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
    'eos8': [
        'Poste grand calcul Onera-ld (ONERA) avec Centos8',
        'ifort', # f77compiler
        'ifort', # f90compiler
        'icc', # Cppcompiler
        ['-DCACHELINE=32', '-DNB_SOCKET=1'], # CppAdditionalOptions
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
    'eosz': [
        'Poste grand calcul eosXXXz (ONERA)',
        'ifort', # f77compiler
        'ifort', # f90compiler
        'icc', # Cppcompiler
        [], # CppAdditionalOptions
        [], # f77AdditionalOptions
        True, # useOMP
        False, # static
        ['/usr/local/hdf5-gnu-1.8.8/include'], # additionalIncludePaths
        [], # additionalLibs
        ['/usr/local/hdf5-gnu-1.8.8/lib'], # additionalLibPaths
        False, # useCuda
        [] # NvccAdditionalOptions
    ],
    ###############################################################################
    'eos': [
        'Onera-eos (legacy-doit etre apres eosZ)',
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
    ],
    ###############################################################################
    'ld': [
        'Poste grand calcul Onera-ld (ONERA)',
        'gfortran', # f77compiler
        'gfortran', # f90compiler
        'gcc', # Cppcompiler
        ['-DCACHELINE=64', '-DNB_SOCKET=1'], # CppAdditionalOptions
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
    'mangrove': [
        'Machine avec acces GPU (ONERA)',
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
    ],
    ###############################################################################
    'moloch': [
        'Machine dev Cedre (ONERA)',
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
    ],
    ###############################################################################
    'btmclx2': [
        'Cluster Turbomeca (Safran)',
        'gfortran', # f77compiler
        'gfortran', # f90compiler
        'gcc', # Cppcompiler
        [], # CppAdditionalOptions
        [], # f77AdditionalOptions
        True, # useOMP
        False, # static
        [], # additionalIncludePaths
        [], # additionalLibs
        ['/usr/lib/gcc/x86_64-redhat-linux/4.1.2'], # additionalLibPaths
        False, # useCuda
        [] # NvccAdditionalOptions
    ],
    ###############################################################################
    'visio': [
        'Machine de post gfx (ONERA)',
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
    'visung': [
        'Machine de post gfx (ONERA)',
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
    ],
    ##############################################################################
    'visung_el8': [
        'Machine de post gfx (ONERA)',
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
    'elmer': [
        'Machine de gros post gfx (ONERA)',
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
    'fulvio': [
        'Machine post gfx legacy (ONERA)',
        'ifort', # f77compiler
        'ifort', # f90compiler
        'icc', # Cppcompiler
        [], # CppAdditionalOptions
        [], # f77AdditionalOptions
        True, # useOMP
        False, # static
        ['/usr/local/hdf5-intel-1.8.8/include'], # additionalIncludePaths
        [], # additionalLibs
        ['/usr/local/hdf5-intel-1.8.8/lib'], # additionalLibPaths
        False, # useCuda
        [] # NvccAdditionalOptions
    ],
    ###############################################################################
    'cobalt': [
        'CCRT machine Cobalt',
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
    'irene': [
        'TGCC machine Irene-skylake-partition',
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
    'topaze': [
        'CCRT machine Topaze-milan-partition',
        'ifort', # f77compiler
        'ifort', # f90compiler
        'icc', # Cppcompiler
        ['-DCACHELINE=32', '-DNB_SOCKET=2', '-DCORE_PER_SOCK=64', '-axCORE-AVX2', '-mavx2'], # CppAdditionalOptions
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
    'ouessant': [
        'Machine IDRIS IBM  POWER + NVIDIA P100)',
        'pgf90', # f77compiler
        'pgf90', # f90compiler
        'pgcc', # Cppcompiler
        [], # CppAdditionalOptions
        [], # f77AdditionalOptions
        True, # useOMP
        False, # static
        [], # additionalIncludePaths
        ["pgf90", "pgf902", "pgc", "pgmath", "pgkomp", 'omp', "pgf90_rpm1", "rt", "pgf90rtl"], # additionalLibs
        ["/opt/pgi/linuxpower/18.4/lib/"], # additionalLibPaths
        False, # useCuda
        [] # NvccAdditionalOptions
    ],
    ###############################################################################
    'jean_zay': [
        'Machine IDRIS intel + NVIDIA V100)',
        'nvfortran', # f77compiler
        'nvfortran', # f90compiler
        'nvc', # Cppcompiler
        [], # CppAdditionalOptions
        [], # f77AdditionalOptions
        True, # useOMP
        False, # static
        [], # additionalIncludePaths
        [], # additionalLibs
        [], # additionalLibPaths
        True, # useCuda
        [] # NvccAdditionalOptions
    ],
    ###############################################################################
    'adastra': [
        'Machine CINES Cray',
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
    'adastra_gpu': [
        'Machine CINES Cray',
        'ftn', # f77compiler
        'ftn', # f90compiler
        'cc', # Cppcompiler
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
    'ubuntu': [
        'Linux ubuntu 24.04',
        'gfortran', # f77compiler
        'gfortran', # f90compiler
        'gcc', # Cppcompiler
        [], # CppAdditionalOptions
        [], # f77AdditionalOptions
        True, # useOMP
        False, # static
        ['/usr/include', '/usr/include/hdf5/openmpi', '/usr/lib/x86_64-linux-gnu/openmpi/include'], # additionalIncludePaths
        ['gfortran', 'gomp'], # additionalLibs
        ['/usr/lib/x86_64-linux-gnu/hdf5/openmpi', '/usr/lib/x86_64-linux-gnu'], # additionalLibPaths
        False, # useCuda
        [] # NvccAdditionalOptions
    ],
    ###############################################################################
    'azure': [
        'Linux Centos7 - Github Actions',
        'gfortran', # f77compiler
        'gfortran', # f90compiler
        'gcc', # Cppcompiler
        [], # CppAdditionalOptions
        [], # f77AdditionalOptions
        True, # useOMP
        False, # static
        ['/usr/include', '/usr/include/hdf5/openmpi', '/usr/lib/x86_64-linux-gnu/openmpi/include'], # additionalIncludePaths
        ['gfortran', 'gomp'], # additionalLibs
        ['/usr/lib/x86_64-linux-gnu/hdf5/openmpi', '/usr/lib/x86_64-linux-gnu'], # additionalLibPaths
        False, # useCuda
        [] # NvccAdditionalOptions
    ],
    ###############################################################################
    'linux': [
        'Linux (generic)',
        'gfortran', # f77compiler
        'gfortran', # f90compiler
        'gcc', # Cppcompiler
        [], # CppAdditionalOptions
        [], # f77AdditionalOptions
        True, # useOMP
        False, # static
        ['/usr/include', '/usr/include/hdf5/openmpi', '/usr/lib/x86_64-linux-gnu/openmpi/include'], # additionalIncludePaths
        ['gfortran', 'gomp'], # additionalLibs
        ['/usr/lib/x86_64-linux-gnu/hdf5/openmpi', '/usr/lib/x86_64-linux-gnu'], # additionalLibPaths
        False, # useCuda
        [] # NvccAdditionalOptions
    ],
    ###############################################################################
    'pdev': [
        'Machine Airbus (Airbus)',
        'ifort', # f77compiler
        'ifort', # f90compiler
        'icc', # Cppcompiler
        [], # CppAdditionalOptions
        [], # f77AdditionalOptions
        True, # useOMP
        False, # static
        ['/opt/soft/cdtng/tools/portage/1.9/usr/include', '/opt/hpmpi/include'], # additionalIncludePaths
        [], # additionalLibs
        ['/opt/soft/cdtng/tools/portage/1.9/usr/lib', '/opt/hpmpi/lib', '/opt/soft/cdtng/tools/intelcompiler/11.0/lib/intel64'], # additionalLibPaths
        False, # useCuda
        [] # NvccAdditionalOptions
    ],
    ###############################################################################
    'sator_coda2': [
        'Cluster de calcul Sator Saphire (ONERA)',
        'gfortran', # f77compiler
        'gfortran', # f90compiler
        'gcc', # Cppcompiler
        ['-DCACHELINE=64', '-DNB_SOCKET=2', '-DCORE_PER_SOCK=48', '-DSIMD=AVX2P512'], # CppAdditionalOptions
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
    'sator': [
        'Cluster de calcul Sator Broadwell (ONERA)',
        'ifort', # f77compiler
        'ifort', # f90compiler
        'icc', # Cppcompiler
        ['-DCACHELINE=32', '-DNB_SOCKET=2', '-DCORE_PER_SOCK=14', '-Dvtune', '-DSIMD=AVX2'], # CppAdditionalOptions
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
    'sator_sky': [
        'Cluster de calcul Sator Skylake (ONERA)',
        'ifort', # f77compiler
        'ifort', # f90compiler
        'icc', # Cppcompiler
        ['-DCACHELINE=64', '-DNB_SOCKET=2', '-DCORE_PER_SOCK=22', '-Dvtune', '-DSIMD=AVX2P512'], # CppAdditionalOptions
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
    'sator_cas': [
        'Cluster de calcul Sator Cascadelake (ONERA)',
        'ifort', # f77compiler
        'ifort', # f90compiler
        'icc', # Cppcompiler
        ['-DCACHELINE=64', '-DNB_SOCKET=2', '-DCORE_PER_SOCK=24', '-Dvtune', '-DSIMD=AVX2P512'], # CppAdditionalOptions
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
    'sator_sph': [
        'Cluster de calcul Sator Saphire (ONERA)',
        'ifort', # f77compiler
        'ifort', # f90compiler
        'icc', # Cppcompiler
        ['-DCACHELINE=64', '-DNB_SOCKET=2', '-DCORE_PER_SOCK=48', '-DSIMD=AVX2P512'], # CppAdditionalOptions
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
    'sator_gcc': [
        'Cluster de calcul Sator Saphire (ONERA)',
        'gfortran', # f77compiler
        'gfortran', # f90compiler
        'gcc', # Cppcompiler
        ['-DCACHELINE=64', '-DNB_SOCKET=2', '-DCORE_PER_SOCK=48', '-DSIMD=AVX2P512'], # CppAdditionalOptions
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
    'juno_gpu': [
        'GPU A30 onera',
        'nvfortran', # f77compiler
        'nvfortran', # f90compiler
        'nvc', # Cppcompiler
        [], # CppAdditionalOptions
        [], # f77AdditionalOptions
        True, # useOMP
        False, # static
        [], # additionalIncludePaths
        [], # additionalLibs
        [], # additionalLibPaths
        True, # useCuda
        [] # NvccAdditionalOptions
    ],
    ###############################################################################
    'juno_gcc': [
        'Machine dev Juno rocky8 (ONERA)',
        'gfortran', # f77compiler
        'gfortran', # f90compiler
        'gcc', # Cppcompiler
        ['-DCACHELINE=64','-DNB_SOCKET=2','-DCORE_PER_SOCK=48','-DSIMD=AVX512'], # CppAdditionalOptions
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
    'juno_coda': [
        'Machine dev Juno rocky8 (ONERA) (env. coda)',
        'gfortran', # f77compiler
        'gfortran', # f90compiler
        'gcc', # Cppcompiler
        ['-DCACHELINE=64','-DNB_SOCKET=2','-DCORE_PER_SOCK=48','-DSIMD=AVX512'], # CppAdditionalOptions
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
    'juno': [
        'Machine dev Juno rocky8 (ONERA)',
        'ifort', # f77compiler
        'ifort', # f90compiler
        'icc', # Cppcompiler
        ['-DCACHELINE=64','-DNB_SOCKET=2','-DCORE_PER_SOCK=48','-DSIMD=AVX512'], # CppAdditionalOptions
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
    'spiro': [
        'Machine dev Spiro (ONERA)',
        'ifort', # f77compiler
        'ifort', # f90compiler
        'icc', # Cppcompiler
        ['-DCACHELINE=32','-DNB_SOCKET=2','-DCORE_PER_SOCK=12','-DSIMD=AVX2'], # CppAdditionalOptions
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
    'chi85bi': [
        'Cluster EDF (ONERA)',
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
    'Raspail': [
        'Machine DTIM (ONERA)',
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
        ['/usr/lib/gcc/x86_64-linux-gnu/7', '/usr/lib/x86_64-linux-gnu/'],  # additionalLibPaths
        False, # useCuda
        [] # NvccAdditionalOptions
    ],
    ###############################################################################
    'curie': [
        'Cluster Curie',
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
    ###############################################################################
    'localhost.localdomain': [
        'Unknown',
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
    'dumbo': [
        'Grosse machine de post-traitement (ONERA)',
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
    'default': [
        'Default',
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
}
