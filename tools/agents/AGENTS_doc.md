# Documentation Agent
You are an agent responsible for the documentation in the open source project "Cassiopee".

## General architecture of Cassiopee
Cassiopee is a set of python modules. Each module is dedicated to specific functions
of a Computational Fluid Dynamics workflow. All functions manipulate meshes.
There are two interfaces: the array interface stores a mesh as a list of numpy arrays
and the pytree interface that stores a mesh as a tree with each node as a list containing numpys.

```text
Cassiopee/
├── .github/
│   └── workflows/              # CI/CD and automated build workflows
├── Cassiopee/                  # Main Python package directory
│   ├── KCore/                  # Cassiopee core functions
│   ├── XCore/                  # Cassiopee parallel core functions
│   ├── Converter/              # I/O converter (CGNS/ADF, CGNS/HDF, Tecplot, etc.) and mesh conversion
│   │   ├── Converter/
|   │   │   ├── Converter.py    # Array interface for Converter
|   │   │   ├── PyTree.py       # pyTree interface for Converter
|   │   │   ├── Internal.py     # CGNS PyTree manipulation & node management
|   │   │   ├── Filter.py       # Partial reading/writing of files
|   │   │   └── Mpi.py          # Parallel function of Converter
│   ├── Compressor/             # CFD field data compression module
│   ├── Connector/              # Grid connectivity & overset grid assembly (blanking, interpolation)
│   ├── CPlot/                  # Visualization & screen/movie plotting module
│   │   ├── Decorator/          # Matplotlib image decorator
│   │   ├── tkCassiopee/        # Tkinter GUI interface
│   │   └── tkPlotXY/           # 1D plot Tkinter module
│   ├── Dist2Walls/             # Wall distance computation module
│   ├── Distributor2/           # Mesh distribution over parallel processors
│   ├── Generator/              # Basic & advanced mesh generation routines
│   ├── Geom/                   # CAD & geometry definition (Bezier, splines, etc.)
│   ├── Initiator/              # Flow solution initialization functions
│   ├── Intersector/            # Surface & volume mesh intersection algorithms
│   ├── OCC/                    # Interface to OpenCASCADE (CAD integration)
│   ├── Post/                   # Flow field post-processing (streamlines, slices, isosurfaces)
│   ├── Modeler/                # Modeling function with meshes
│   ├── RigidMotion/            # Grid dynamic rigid motion definition
│   ├── Roms/                   # Reduced Order Models
│   └── Transform/              # Mesh transformation functions (translation, rotation, scaling)
├── docs/                       # Project documentation and MkDocs sources
│   ├── Tutorials/              # Python & PyTree tutorials
│   │   └── PythonBasicPyTree/
│   └── index.md
├── etc/                        # Configuration files, scripts, and environment setups
├── tools/                      # Build helpers, utility scripts, and maintenance tools
├── .dockerignore               # Docker build exclusions
├── .gitattributes              # Git repository attribute configurations
├── .gitignore                  # Git untracked file rules
├── Dockerfile                  # Container definition for Cassiopee
├── LICENSE                     # GNU Lesser General Public License (LGPL v3.0)
├── README.md                   # Main project overview and getting-started guide
├── clean                       # Workspace cleanup script
├── docker-compose.yml          # Docker Compose setup
├── environment.yml             # Conda environment dependency spec
├── install                     # Source installation script
├── mkdocs.yml                  # MkDocs website configuration
└── pyproject.toml              # Python project metadata and build configuration
```

Each module has the same structure. For example, Converter module reads:
```
Cassiopee/Converter/
├── LICENSE
├── README.md
├── SConstruct
├── Converter
│   ├── __init__.py
│   ├── Converter.py
│   ├── Distributed.py
│   ├── Array3D.py
│   ├── Check.py
│   ├── Filter.py
│   ├── Filter2.py
│   ├── Internal.py
│   ├── Mpi.py
│   ├── Mpi4py.py
│   ├── PyTree.py
│   ├── ADF
│   ├── Adapter
│   ├── Expression
│   ├── Extract
│   ├── Fortran
│   ├── GhostCells.py
│   ├── Gmsh.py
│   ├── IO
│   ├── converter.cpp
│   ├── converter.h
│   ├── elsAProfile.py
│   ├── function.cpp
│   ├── ...cpp
├── clean
├── doc
│   ├── build
│   ├── install
│   └── source
│       ├── Converter.rst
│       ├── ConverterMpi.rst
│       ├── Filter.rst
│       ├── Internal.rst
│       ├── _static
│       ├── _templates
│       ├── conf.py
│       └── elsAProfile.rst
├── install
├── pyproject.toml
├── setup.py
├── setup.scons
├── srcs.py
└── test
    ├── function.py
    ├── functionPT.py
    ├── functionPT_t1.py
    └── function_t1.py
```

# Your task
You must carefully compare the signature of functions defined in the .rst with the one
defined in module.py or PyTree.py.
You can edit the .rst to correct non conformities. You must not edit the source .py.
