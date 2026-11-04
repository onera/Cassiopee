<p align="center">
  <img src="./docs/Images/CassiopeeTitreWhiteBg.png" width="80%" alt="Cassiopee Logo">
</p>

**A Python toolkit for CFD pre-processing, mesh generation, and post-processing.**


[![Documentation](https://img.shields.io/badge/Documentation-online-blue)](https://onera.github.io/Cassiopee/latest/)
[![Latest release](https://img.shields.io/github/v/release/onera/Cassiopee?label=release)](https://github.com/onera/Cassiopee/releases)
[![License: LGPL-3.0](https://img.shields.io/badge/License-LGPL--3.0-blue.svg)](LICENSE)
[![GitHub stars](https://img.shields.io/github/stars/onera/Cassiopee?style=social)](https://github.com/onera/Cassiopee/stargazers)


Cassiopee is a collection of 19 Python modules developed by **ONERA, the French Aerospace Lab**, for computational fluid dynamics (CFD) workflows.

Its modular API provides building blocks for geometry processing, mesh generation and transformation, data conversion, and solution post-processing. Cassiopee is solver-agnostic and supports workflows involving structured and unstructured meshes, including advanced overset-grid and immersed-boundary methods.


**Explore the documentation, run an example, ask for help, or install Cassiopee using the links below.**

- [Documentation](https://onera.github.io/Cassiopee/latest/)
- [Tutorials](https://github.com/onera/Cassiopee/blob/dev/docs/Tutorials/Tutorials.md) (work in progress)
- [Release notes](https://onera.github.io/Cassiopee/latest/ReleaseNotes.html)
- [Report a bug or ask a question](https://github.com/onera/Cassiopee/issues)

<br>

## What can you do with Cassiopee?

Cassiopee provides specialized modules that can be combined into scripts and application-specific workflows. Most functionalities can run efficiently in parallel, on both shared-memory and distributed-memory systems.

| Module | Main capabilities |
|---|---|
| **Converter** | Read, write, and manipulate CFD data in different formats, including CGNS and Tecplot. |
| **Geom** | Create and manipulate geometric entities. |
| **Generator** | Generate structured and unstructured meshes. |
| **Transform** | Apply transformations to meshes. |
| **Intersector** | Perform mesh intersection operations. |
| **Connector** | Build mesh connectivity, including overset-grid connectivity and immersed-boundary data via a dedicated submodule. |
| **Initiator** | Initialize solution fields. |
| **Post** | Extract and post-process CFD results, including slices, streamlines, and isosurfaces. |
| **Dist2Walls** | Compute distances to walls and bodies. |
| **Distributor2** | Distribute mesh blocks across processors. |
| **CPlot** | Visualize meshes and solution fields via the `cassiopee` and `cplot` graphical user interfaces (GUIs). |
| **OCC** | Interface with Open CASCADE for CAD-related operations. |

Additional modules and submodules provide capabilities such as rigid-body motion, field compression, and parallel data management. See the [complete documentation](https://onera.github.io/Cassiopee/latest/) for more details.

<br>

## Designed for Python workflows

Cassiopee offers two complementary data interfaces:

- **Arrays:** a lightweight mesh and data representation based on numpy arrays.
- **CGNS/Python trees (pyTrees):** a hierarchical representation of numerical data built from lists and numpy arrays, following the CFD General Notation System (CGNS) standard.

Both interfaces can be used directly from Python scripts to automate repetitive tasks and build reusable CFD workflows.

### A first example

The following script generates a Cartesian mesh and writes it to a CGNS file:

```python
import Generator.PyTree as G
import Converter.PyTree as C

# Generate a Cartesian mesh
zone = G.cart((0., 0., 0.), (0.1, 0.1, 0.1), (11, 11, 11))

# Create a CGNS tree and write it to disk
tree = C.newPyTree(["Base", zone])
C.convertPyTree2File(tree, "mesh.cgns")
```

For more examples, explore the [tutorials](https://github.com/onera/Cassiopee/blob/dev/docs/Tutorials/Tutorials.md) and the tests provided for each Cassiopee module.

<br>

## Installation

Choose the installation method that best fits your environment.

### Install from binaries

#### Linux — pip

```bash
python -m pip install cassiopee  # int32
python -m pip install cassiopee-i8  # int64
```

Binary availability depends on the Python version and target platform. See the installation documentation if no compatible package is available for your setup.

#### Windows

See the available [Windows releases](https://github.com/onera/Cassiopee/releases).

#### Docker

For a container-based installation, follow the [Docker instructions](https://github.com/onera/Cassiopee/blob/main/docs/developers/Docker/UsingDocker.md).

### Build from source

Building from source is supported on the following platforms:

- [Windows — MSYS2](https://github.com/onera/Cassiopee/blob/dev/docs/developers/Install/msys2.md)
- [Linux — Ubuntu](https://github.com/onera/Cassiopee/blob/dev/docs/developers/Install/ubuntu.md)
- Apple — macOS (coming soon)

Building from source may require a suitable compiler toolchain and external libraries. Follow the platform-specific instructions for the required dependencies and environment setup.

<br>

## Executables

| Application | Description |
|---|---|
| `cassiopee` | Main graphical user interface (GUI) for pyTree visualization and pre-/post-processing. |
| `cplot` | Cassiopee's lightweight display application, supporting a wide range of mesh and visualization formats. |
| `validCassiopee` | Cassiopee's continuous integration (CI) tool: runs regression tests and CFD validation cases, comparing results against local/global references and tracking test status and execution times. |
| `kpython` | Cassiopee's Python executable wrapper: facilitates the launch of Python scripts sequentially or in parallel via MPI and OpenMP. Additional options include profiling, sanitizer, and GNU Debugger (GDB) support. |
| `kdiff` | Cassiopee's own diff tool: compares two files across a wide range of file formats. |

<br>

## License

Cassiopee is distributed under the **GNU Lesser General Public License v3.0 (LGPL-3.0)**. See [LICENSE](LICENSE) for the full license text.

<br>

<p align="center">
  <img src="./docs/Images/DroneAD.png" width="80%" alt="CFD flow visualization">
</p>