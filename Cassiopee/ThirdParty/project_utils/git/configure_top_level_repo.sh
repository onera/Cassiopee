#/bin/bash

source git/submodule_utils.sh

git_root_dir="$(git rev-parse --show-superproject-working-tree)"

echo $git_root_dir
cd "${git_root_dir}"
git submodule update --init
git_config_submodules
