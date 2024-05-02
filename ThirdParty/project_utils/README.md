## Introduction ##
`project_utils` offers some git submodule (alias and commands) and cmake utilities (macros functions and Find.cmake scripts). It also documents the intended workflow to use these utilities.

## Dependencies ##
`project_utils` is a cmake 3.12 library.

## Acknowledgements ##
Some of the Find.cmake scripts have been copy-pasted from the internet, then modified. A first-line comment gives the source it was copied from.

## License ##
`project_utils` is available under the MPL-2.0 license (https://mozilla.org/MPL/2.0/).

## TODO ##
TOC:
git and submodules
cmake

build doc with cmake, sphinx, doxygen, breathe
https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/


1. rendre les dépendances optionnelles (et possiblement configurables).
 Pour cela, il va falloir étudier la faisabilité d'une approche propre en sélectionnant l' (ou les) outil(s) le(s) plus adapté(s) à la tâche : git / cmake / spack.
 Si avec cette nouvelle gestion des dépendances, on pouvait éviter le travail manuel de tirer les dépendances récursives, ce serait cool aussi.
2. bootstrap project_utils
 c'est probable qu'une solution simple au problème actuel avec la macro check_sub_repo_exists serait de la mettre dans project_utils (la version optionnelle), et de mettre au début de tout les CMakeList.txt le contrôle de project_utils.
 Important : ceci est indépendant de l'outil utilisé pour la gestion des dépendances, car pour le moment, on garde cmake comme outil de build
3. Etude d'une simplification/refactoring/meilleur documentation de project_utils
 Comme on pense faire quelques modification sur l'env de dev, le projet coeur de cet env sera certainement affecté, et c'est une bonne occasion de voir si des amélioration sont possibles. 
