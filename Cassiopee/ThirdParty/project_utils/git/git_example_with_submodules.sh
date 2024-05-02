#!/bin/bash

ROOT_DIR=$(pwd)

# 0. create server repos
mkdir server
(cd server
  mkdir base_lib.git
  (
    cd base_lib.git
    git init --bare
  )
  mkdir app_lib.git
  (
    cd app_lib.git
    git init --bare
  )
  mkdir my_appli.git
  (
    cd my_appli.git
    git init --bare
  )
)

# 1. create local repos
mkdir local 
(cd local 
  git clone ${ROOT_DIR}/server/base_lib.git
  git clone ${ROOT_DIR}/server/app_lib.git
  git clone ${ROOT_DIR}/server/my_appli.git
)

# 2. commits my_appli
(cd local/my_appli
# 2.0. simple commit
  echo pouet > file0.txt
  git add -A
  git commit -m "pouet in file 0"
  git push
# 2.1. branch
  git checkout -b dev_mxpl
  echo "mxpl work" > mxpl.cpp
  git add -A
  git commit -m "begin mxpl"
  echo "more mxpl work" >> mxpl.cpp
  git add -A
  git commit -m "continue mxpl"
# 2.2. come back to master
  git checkout master
  echo "master work" > master_work.cpp
  git add -A
  git commit -m "master work"
  echo "other master work" > other_master_work.cpp
  git add -A
  git commit -m "other master work"
# 2.3. end work on branch
  git checkout dev_mxpl
  echo "end mxpl work" >> mxpl.cpp
  git add -A
  git commit -m "super mxpl working"
# 2.4. merge it to trunk
  git checkout master
  git status
  git merge dev_mxpl
  git status
# 2.5. push
  git push
)
  
# 2. commits libs
(cd local/base_lib
  echo "my base algo" > algo.hpp
  git status
  git add -A
  git commit -m "algo"
  git push
)
(cd local/app_lib
  echo "my app algo" > algo.hpp
  git status
  git add -A
  git commit -m "algo"
  git push
)

# 3. use libs as a submodules
(cd local/app_lib
  git submodule add ${ROOT_DIR}/server/base_lib.git/ external/base_lib
  git commit -m "added submodule"
  git push
)
(cd local/my_appli
  git submodule add ${ROOT_DIR}/server/base_lib.git/ external/base_lib
  git submodule add ${ROOT_DIR}/server/app_lib.git/ external/app_lib
  git commit -m "added submodules"
  git push
)

# 4. start fresh
(cd local
  rm -rf base_lib
  rm -rf app_lib
  rm -rf my_appli
  git sclone ${ROOT_DIR}/server/my_appli.git
)

# 5. dev within 3 repositories
(cd local/my_appli
  # 5.0. same branch
  cd external/base_lib
  echo "new base algo" > new_algo.hpp
  git add -A
  git commit -m "new base algo"

  cd ../app_lib
  echo "new app algo" > new_algo.hpp
  git add -A
  git commit -m "new app algo"

  cd ../..
  echo "#include <base_lib/new_algo.hpp>" > use.cpp 
  echo "#include <app_lib/new_algo.hpp>" >> use.cpp 
  git status
  git add -A
  git commit -m "use of new algo"
  git spush

  # 5.1. several branches
  cd external/base_lib
  git tag v1.0
  git push origin v1.0

  cd ../..
  git checkout -b my_dev

  cd external/base_lib
  echo "my_dev base algo" > my_dev_algo.hpp
  rm new_algo.hpp
  git add -A
  git commit -m "my_dev base algo"

  cd ../app_lib
  echo "my_dev algo" > new_algo.hpp
  git add -A
  git commit -m "my_dev app algo"

  cd ../..
  echo "#include <base_lib/my_dev_algo.hpp>" > use.cpp 
  echo "#include <app_lib/new_algo.hpp>" >> use.cpp 
  echo "[using my_dev algos]" >> use.cpp 
  git status
  git add -A
  git commit -m "use of my_dev algo"
  git spush origin my_dev

  # 5.2. come back to master (with coherent submodules working tree)
  git scheckout master
)



