# Template
A template module

This module is just a template module.
You must first adapt it to your module name:
python rename.py <Module>
where <Module> is the name of your module.

In Template, you will find an example of use of Cassiopee array interface
in arrayExample.cpp.
You will also find an example of use of pyTree interface in pyTreeExample.cpp.

This module is only dependant on Cassiopee/KCore module (LGPL3).
Feel free to use it as a canvas.

For support, go to https://github.com/onera/Cassiopee

## installing this module:

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
