# Generator
Generating grids

ONERA.
Licenced under LGPL3.

For support, go to https://github.com/onera/Cassiopee

This distribution contains third party software:
- netgen (Generator/Netgen, see Generator/Netgen/LICENCE)
- tetgen (Generator/Tetgen, see Generator/Tetgen/LICENCE)

## installing Generator:

- before installing, edit etc/installBaseUser.json. Set the correct 
compiler names and libraries.

Then,

For standard install:
- set CASSIOPEE variable to Cassiopee directory, for instance:
  in sh/bash: export CASSIOPEE="/home/john/Cassiopee"
  in csh/tcsh: setenv CASSIOPEE "/home/john/Cassiopee"
- set ELSAPROD variable to your processor type, for instance:
  in sh/bash: export ELSAPROD="x86"
  in csh/tcsh: setenv ELSAPROD "x86"
- type: ./install
- binaries will be installed to $CASSIOPEE/Dist/bin/$ELSAPROD
- follow instructions

For custom install:
- set CASSIOPEE variable to Cassiopee source directory, for instance:
  in sh/bash: export CASSIOPEE="/home/john/Cassiopee"
  in csh/tcsh: setenv CASSIOPEE "/home/john/Cassiopee"
- set ELSAPROD variable to your processor type, for instance:
  in sh/bash: export ELSAPROD="x86"
  in csh/tcsh: setenv ELSAPROD "x86"
- set CASSIOPEE_INSTALL_PREFIX to the target directory for install.
  in sh/bash: export CASSIOPEE_INSTALL_PREFIX="/home/john/<YOUR_DIRECTORY>"
  in csh/tcsh: setenv CASSIOPEE_INSTALL_PREFIX "/home/john/<YOUR_DIRECTORY>"

- type: ./install
- binaries will be installed to <YOUR_DIRECTORY>
- follow instructions

## Using Cassiope modules:
For standard install
- in sh/bash: source $CASSIOPEE/Dist/env_Cassiopee.sh
  in csh/tcsh: source $CASSIOPEE/Dist/env_Cassiopee.sh
- run "python" or "cassiopee" (if you have CPlot installed)

For custom install:
- in sh/bash: source <YOUR_DIRECTORY>/env_Cassiopee.sh
  in csh/tcsh: source <YOUR_DIRECTORY>/env_Cassiopee.csh
- run "python" or "cassiopee" (if you have CPlot installed)
