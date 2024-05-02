## Main commands ##

# more info on submodules when doing a git status
git config --global status.submoduleSummary true

# git diff with more verbose info on the submodule changes
git config --global diff.submodule log # equivalent to on-spot `git diff --submodule=log`

# add submodule
git submodule add <remote of submodule> <sub-directory where the submodule will be put>

# clone repository with submodules
git clone <my_repository/my_project>
cd my_project
git submodule update --init # init: tell git to track submodules, update: ~= for each submodule s, git clone s
git status # to ensure all is OK

# git pull with submodules
git pull # pull the container
git submodule sync # update submodules remote location
git submodule update # pull the submodules


## Warnings ##
# 1. Do NOT use the `recursive` flag in:
git clone
git submodule update --init
# it is not what we want to do with the current workflow, 
# since each repository is supposed to look at its direct submodules NOT recursively,
# because it is supposed to have its nested submodules as direct submodules
# Example:
# if sub-B is submodule of A
# if sub-sub-C is submodule of sub-B
# then sub-sub-C MUST ALSE BE a submodule of A.
# hence we do not want to confuse the build system, or ourself as a developper, into having sub-sub-C checked out more than once


## Troubleshooting ##

# if after `git submodule update`
# with `git status` you see 
#     modified:   external/sub-A (modified content)
# then it means something went wrong with sub-A. Fix it with
cd external/mod-A
git status
# if it shows something like
#     On branch my_branch
#     Your branch is up-to-date with 'origin/my_branch'.
#  
#     Changes to be committed:
#       (use "git reset HEAD <file>..." to unstage)
#     
#             deleted:    file1
#             deleted:    file2
#             ...
# The working tree is not there
# So simply check it out
git checkout my_branch .  # the dot tells git you are not going to another branch, only copying the files to the working tree
