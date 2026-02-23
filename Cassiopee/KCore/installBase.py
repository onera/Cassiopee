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

    'WDAAA728Z': {
        'description': 'Windows win64+msys2 (XJ-ONERA)',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': [
            '-D__SHADERS__',
            '-isystem /d/juvigny/msys64/mingw64/include/python3.8',
            '-isystem /d/juvigny/msys64/mingw64/lib/python3.8/site-packages/numpy/core/include/numpy/'
        ],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [
            'd:/juvigny/msys64/mingw64/include',
            'd:/juvigny/msys64/mingw64/include/OpenBLAS'
        ],
        'additionalLibs': ['gfortran', 'gomp', 'pthread', 'openblas', 'psapi'],
        'additionalLibPaths': [
            'd:/juvigny/msys64/mingw64/lib',
            'd:/juvigny/msys64/mingw64/bin'
        ],
        'useCuda': True,
        'NvccAdditionalOptions': ['-arch=sm_60']
    },

    'Aryen': {
        'description': 'Windows win64+msys2 (CB-Home)',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': ['-Wno-attributes', '-DSIMD=AVX2'],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': ['c:/msys64/mingw64/include'],
        'additionalLibs': ['gomp', 'gfortran'],
        'additionalLibPaths': ['c:/msys64/mingw64/lib'],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'win64': {
        'description': 'Windows win64+msys2',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': ['-Wno-attributes'],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': ['gomp', 'gfortran'],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'WDAAA161Z': {
        'description': 'Windows win64+msys2 (CB-ONERA)',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': ['-Wno-attributes', '-fcommon'],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [
            'c:/Users/benoit/msys64/mingw64/include',
            'c:/Users/benoit/msys64/mingw64/include/openblas'
        ],
        'additionalLibs': ['gomp', 'gfortran'],
        'additionalLibPaths': ['c:/Users/benoit/msys64/mingw64/lib'],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'node6.cluster': {
        'description': 'MacOSX (generic)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': ['/usr/X11/include'],
        'additionalLibs': ['python2.7', 'ifcore'],
        'additionalLibPaths': [
            '/usr/X11/lib',
            '/System/Library/Frameworks/OpenGL.framework/Libraries/'
        ],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'd1log1': {
        'description': 'Cluster HPC4B dev/val (Airbus)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icpc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [
            '/opt/hpmpi/include',
            '/opt/soft/cdtng/tools/portage/1.11/usr/include'
        ],
        'additionalLibs': ['svml', 'irc', 'ifcore', 'iomp5'],
        'additionalLibPaths': [
            '/opt/soft/cdtng/tools/intelcompiler/16.0/compiler/lib/intel64',
            '/opt/soft/cdtng/tools/portage/1.11/usr/lib',
            '/opt/hpmpi/lib/linux_amd64'
        ],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'devt1n007.': {
        'description': 'Cluster HPC5 dev/val tou_b (Airbus)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icpc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [
            '/opt/mpi10/include',
            '/opt/soft/cdtng/tools/portage/1.11/usr/include'
        ],
        'additionalLibs': ['ifcore', 'iomp5', 'svml', 'irc'],
        'additionalLibPaths': [
            '/opt/soft/cdtng/tools/intelcompiler/2018/compilers_and_libraries_2018.5.274/linux/compiler/lib/intel64',
            '/opt/soft/cdtng/tools/portage/1.11/usr/lib',
            '/opt/mpi10/lib'
        ],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'AA': {
        'description': 'Cluster AA',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': ['-DCACHELINE=16'],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [
            '/soft/intel-2017.5/compilers_and_libraries_2017.5.239/linux/mpi/include64',
            '/home/mtaplaf/Environnement_Commun/Elsa/LIB_EXT/HDF5/hdf5-1.8.18-par/include'
        ],
        'additionalLibs': ['mpi'],
        'additionalLibPaths': [
            '/soft/intel-2017.5/compilers_and_libraries_2017.5.239/linux/mpi/lib64',
            '/home/mtaplaf/Environnement_Commun/Elsa/LIB_EXT/HDF5/hdf5-1.8.18-par/lib'
        ],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'AA2': {
        'description': 'Machine AA2',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': ['-DCACHELINE=16'],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [
            '/soft/intel-2017.5/compilers_and_libraries_2017.5.239/linux/mpi/include64',
            '/home/mtaplaf/Environnement_Commun/Elsa/LIB_EXT/HDF5/hdf5-1.8.18-par/include'
        ],
        'additionalLibs': ['mpi'],
        'additionalLibPaths': [
            '/soft/intel-2017.5/compilers_and_libraries_2017.5.239/linux/mpi/lib64',
            '/home/mtaplaf/Environnement_Commun/Elsa/LIB_EXT/HDF5/hdf5-1.8.18-par/lib'
        ],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'wfrontend1': {
        'description': 'Cluster Kairos (Safran)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': False,
        'static': False,
        'additionalIncludePaths': [
            '/appl/APPLI_SNECMA/HDF5/oper/1.8.11/include'
        ],
        'additionalLibs': ['ifcore', 'svml', 'irc'],
        'additionalLibPaths': [
            '/opt/intel/composer_xe_2013_sp1.0.080/lib/intel64',
            '/appl/APPLI_SNECMA/HDF5/oper/1.8.11/lib'
        ],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'celeste': {
        'description': 'Grosse machine de post-traitement (ONERA)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': ['Xxf86vm'],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'linux64': {
        'description': 'Production linux64 (generic)',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': ['/stck1/benoit/include'],
        'additionalLibs': [],
        'additionalLibPaths': ['/stck1/benoit/lib'],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'ld-clang': {
        'description': 'Poste grand calcul Onera-ld (ONERA)',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'clang',
        'CppAdditionalOptions': ['-DCACHELINE=32', '-DNB_SOCKET=1'],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'ld_spack2': {
        'description': 'Poste grand calcul Onera-ld (ONERA)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': ['-DCACHELINE=32', '-DNB_SOCKET=1'],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'eos8': {
        'description': 'Poste grand calcul Onera-ld (ONERA) avec Centos8',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': ['-DCACHELINE=32', '-DNB_SOCKET=1'],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'eosz': {
        'description': 'Poste grand calcul eosXXXz (ONERA)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': ['/usr/local/hdf5-gnu-1.8.8/include'],
        'additionalLibs': [],
        'additionalLibPaths': ['/usr/local/hdf5-gnu-1.8.8/lib'],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'eos': {
        'description': 'Onera-eos (legacy-doit etre apres eosZ)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'ld': {
        'description': 'Poste grand calcul Onera-ld (ONERA)',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': ['-DCACHELINE=64', '-DNB_SOCKET=1'],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'ld_coda': {
        'description': 'Poste grand calcul Onera-ld (ONERA)',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': ['-DCACHELINE=64', '-DNB_SOCKET=1'],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'mangrove': {
        'description': 'Machine avec acces GPU (ONERA)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'moloch': {
        'description': 'Machine dev Cedre (ONERA)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'btmclx2': {
        'description': 'Cluster Turbomeca (Safran)',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': ['/usr/lib/gcc/x86_64-redhat-linux/4.1.2'],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'visio': {
        'description': 'Machine de post gfx (ONERA)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': ['Xxf86vm'],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'visung': {
        'description': 'Machine de post gfx (ONERA)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'visung_el8': {
        'description': 'Machine de post gfx (ONERA)',
        'f77compiler': 'ifx',
        'f90compiler': 'ifx',
        'Cppcompiler': 'icx',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': ['Xxf86vm'],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'elmer': {
        'description': 'Machine de gros post gfx (ONERA)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': ['Xxf86vm'],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'fulvio': {
        'description': 'Machine post gfx legacy (ONERA)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': ['/usr/local/hdf5-intel-1.8.8/include'],
        'additionalLibs': [],
        'additionalLibPaths': ['/usr/local/hdf5-intel-1.8.8/lib'],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'cobalt': {
        'description': 'CCRT machine Cobalt',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': ['-DCACHELINE=32'],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'irene': {
        'description': 'TGCC machine Irene-skylake-partition',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': ['-DCACHELINE=32'],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'topaze': {
        'description': 'CCRT machine Topaze-milan-partition',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [
            '-DCACHELINE=32',
            '-DNB_SOCKET=2',
            '-DCORE_PER_SOCK=64',
            '-axCORE-AVX2',
            '-mavx2'
        ],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'ouessant': {
        'description': 'Machine IDRIS IBM  POWER + NVIDIA P100)',
        'f77compiler': 'pgf90',
        'f90compiler': 'pgf90',
        'Cppcompiler': 'pgcc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [
            'pgf90', 'pgf902', 'pgc', 'pgmath',
            'pgkomp', 'omp', 'pgf90_rpm1', 'rt', 'pgf90rtl'
        ],
        'additionalLibPaths': ['/opt/pgi/linuxpower/18.4/lib/'],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'jean_zay': {
        'description': 'Machine IDRIS intel + NVIDIA V100)',
        'f77compiler': 'nvfortran',
        'f90compiler': 'nvfortran',
        'Cppcompiler': 'nvc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': True,
        'NvccAdditionalOptions': []
    },

    'adastra': {
        'description': 'Machine CINES Cray',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'adastra_gpu': {
        'description': 'Machine CINES Cray',
        'f77compiler': 'ftn',
        'f90compiler': 'ftn',
        'Cppcompiler': 'cc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'ubuntu': {
        'description': 'Linux ubuntu 24.04',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [
            '/usr/include',
            '/usr/include/hdf5/openmpi',
            '/usr/lib/x86_64-linux-gnu/openmpi/include'
        ],
        'additionalLibs': ['gfortran', 'gomp'],
        'additionalLibPaths': [
            '/usr/lib/x86_64-linux-gnu/hdf5/openmpi',
            '/usr/lib/x86_64-linux-gnu'
        ],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'azure': {
        'description': 'Linux Centos7 - Github Actions',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [
            '/usr/include',
            '/usr/include/openmpi-x86_64'
        ],
        'additionalLibs': ['gfortran', 'gomp'],
        'additionalLibPaths': [
            '/usr/lib64',
            '/usr/lib64/openmpi/lib'
        ],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'linux': {
        'description': 'Linux (generic)',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [
            '/usr/include',
            '/usr/include/hdf5/openmpi',
            '/usr/lib/x86_64-linux-gnu/openmpi/include'
        ],
        'additionalLibs': ['gfortran', 'gomp'],
        'additionalLibPaths': [
            '/usr/lib/x86_64-linux-gnu/hdf5/openmpi',
            '/usr/lib/x86_64-linux-gnu'
        ],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'pdev': {
        'description': 'Machine Airbus (Airbus)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [
            '/opt/soft/cdtng/tools/portage/1.9/usr/include',
            '/opt/hpmpi/include'
        ],
        'additionalLibs': [],
        'additionalLibPaths': [
            '/opt/soft/cdtng/tools/portage/1.9/usr/lib',
            '/opt/hpmpi/lib',
            '/opt/soft/cdtng/tools/intelcompiler/11.0/lib/intel64'
        ],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'sator_coda': {
        'description': 'Cluster de calcul Sator Saphire (ONERA)',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': [
            '-DCACHELINE=64',
            '-DNB_SOCKET=2',
            '-DCORE_PER_SOCK=48',
            '-DSIMD=AVX2P512'
        ],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': ['/usr/include/openblas/'],
        'additionalLibs': ['openblas'],
        'additionalLibPaths': ['/usr/lib64/'],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'sator_coda_2511': {
        'description': 'Cluster de calcul Sator Saphire (ONERA)',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': [
            '-DCACHELINE=64',
            '-DNB_SOCKET=2',
            '-DCORE_PER_SOCK=48',
            '-DSIMD=AVX2P512'
        ],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': ['/usr/include/openblas/'],
        'additionalLibs': ['openblas'],
        'additionalLibPaths': ['/usr/lib64/'],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'sator': {
        'description': 'Cluster de calcul Sator Broadwell (ONERA)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [
            '-DCACHELINE=32',
            '-DNB_SOCKET=2',
            '-DCORE_PER_SOCK=14',
            '-Dvtune',
            '-DSIMD=AVX2'
        ],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'sator_sky': {
        'description': 'Cluster de calcul Sator Skylake (ONERA)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [
            '-DCACHELINE=64',
            '-DNB_SOCKET=2',
            '-DCORE_PER_SOCK=22',
            '-Dvtune',
            '-DSIMD=AVX2P512'
        ],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'sator_cas': {
        'description': 'Cluster de calcul Sator Cascadelake (ONERA)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [
            '-DCACHELINE=64',
            '-DNB_SOCKET=2',
            '-DCORE_PER_SOCK=24',
            '-Dvtune',
            '-DSIMD=AVX2P512'
        ],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'sator_sph': {
        'description': 'Cluster de calcul Sator Saphire (ONERA)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [
            '-DCACHELINE=64',
            '-DNB_SOCKET=2',
            '-DCORE_PER_SOCK=48',
            '-DSIMD=AVX2P512'
        ],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'sator_gcc': {
        'description': 'Cluster de calcul Sator Saphire (ONERA)',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': [
            '-DCACHELINE=64',
            '-DNB_SOCKET=2',
            '-DCORE_PER_SOCK=48',
            '-DSIMD=AVX2P512'
        ],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'juno_gpu': {
        'description': 'GPU A30 onera',
        'f77compiler': 'nvfortran',
        'f90compiler': 'nvfortran',
        'Cppcompiler': 'nvc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': True,
        'NvccAdditionalOptions': []
    },

    'juno_gcc': {
        'description': 'Machine dev Juno rocky8 (ONERA)',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': [
            '-DCACHELINE=64',
            '-DNB_SOCKET=2',
            '-DCORE_PER_SOCK=48',
            '-DSIMD=AVX512'
        ],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'juno_coda': {
        'description': 'Machine dev Juno rocky8 (ONERA) (env. coda)',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': [
            '-DCACHELINE=64',
            '-DNB_SOCKET=2',
            '-DCORE_PER_SOCK=48',
            '-DSIMD=AVX512'
        ],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'juno': {
        'description': 'Machine dev Juno rocky8 (ONERA)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [
            '-DCACHELINE=64',
            '-DNB_SOCKET=2',
            '-DCORE_PER_SOCK=48',
            '-DSIMD=AVX512'
        ],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'spiro': {
        'description': 'Machine dev Spiro (ONERA)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [
            '-DCACHELINE=32',
            '-DNB_SOCKET=2',
            '-DCORE_PER_SOCK=12',
            '-DSIMD=AVX2'
        ],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'chi85bi': {
        'description': 'Cluster EDF (ONERA)',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'Raspail': {
        'description': 'Machine DTIM (ONERA)',
        'f77compiler': 'gfortran-7',
        'f90compiler': 'gfortran-7',
        'Cppcompiler': 'clang++-5.0',
        'CppAdditionalOptions': [
            '-pedantic',
            '-march=native',
            '-Wno-variadic-macros',
            '-Wno-long-long',
            '-g'
        ],
        'f77AdditionalOptions': [
            '-march=native',
            '-fdefault-real-8',
            '-fdefault-double-8'
        ],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': ['/usr/include/hdf5/serial/'],
        'additionalLibs': [],
        'additionalLibPaths': [
            '/usr/lib/gcc/x86_64-linux-gnu/7',
            '/usr/lib/x86_64-linux-gnu/'
        ],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'curie': {
        'description': 'Cluster Curie',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': ['/usr/local/hdf5-1.8.8/include'],
        'additionalLibs': [],
        'additionalLibPaths': ['/usr/local/hdf5-1.8.8/lib'],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'localhost.localdomain': {
        'description': 'Unknown',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'dumbo': {
        'description': 'Grosse machine de post-traitement (ONERA)',
        'f77compiler': 'ifort',
        'f90compiler': 'ifort',
        'Cppcompiler': 'icc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': ['Xxf86vm'],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    },

    'default': {
        'description': 'Default',
        'f77compiler': 'gfortran',
        'f90compiler': 'gfortran',
        'Cppcompiler': 'gcc',
        'CppAdditionalOptions': [],
        'f77AdditionalOptions': [],
        'useOMP': True,
        'static': False,
        'additionalIncludePaths': [],
        'additionalLibs': [],
        'additionalLibPaths': [],
        'useCuda': False,
        'NvccAdditionalOptions': []
    }
}
