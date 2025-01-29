# Contributing to the Cassiopee Github repository

## TLDR

To contribute, you need to create a `Pull Request` from your remote fork webpage: https://github.com/github_username/Cassiopee/. Click on the _Contribute_ button and follow the instructions (give a title to your PR, describe the modifications and review your changes)

<br> </br>


## Forking the upstream repository

- To create your own fork of Cassiopee, please read this [section](https://github.com/onera/Cassiopee/blob/main/docs/developers/Git/cloning-cassiopee.md#clone-or-fork-cassiopee)  

<br> </br>

## Synchronising your fork from the dev branch of the upstream repository

- From your remote fork webpage (https://github.com/github_username/Cassiopee/), click on the _Sync fork_ button

- From a terminal window, follow the instructions given at the end of this [paragraph](https://github.com/onera/Cassiopee/blob/main/docs/developers/Git/cloning-cassiopee.md#clone-or-fork-cassiopee)

<br> </br>

## Stashing

Pulling updates will fail if you have local changes. You can first check the state of your local/working directories with

```sh
git status
```

and if anything is listed, please consider stashing it (ie, saving it). Below is a cheat sheet on stashing

```sh
git stash list  # list all stashes in antichronological order of insertion in the stashing stack

git stash push -m "stashing message to help you remember what this stash is about"  # save your local changes in your local stashing stack

git stash apply stash@{stash_id}  # reinstate the modifications saved in `stash_id` in your working directory. stash_id 0 is the most recent.

git stash drop stash@{stash_id}  # delete stash stash_id - other stashes may be reindexed

git stash pop stash@{stash_id}  # combine git stash apply and git stash drop in one command
```

<br> </br>

## Synchronising your local and working directories from your remote fork

To synchronise your local directory only:
```sh
git fetch origin dev
```

To synchronise your working directory from your local directory:
```sh
git merge
```

To synchronise both directories at once:
```sh
git pull origin dev
```

<br> </br>

## Synchronising your local and working directories, and your remote fork, from the upstream repo

```sh
git pull upstream dev
git push origin dev
```

<br> </br>

## Committing and pushing to your remote fork

1. Committing will commit your local modifications to your local directory from your working directory. 
First inspect the files that are modified:
```sh
git status
```

To see the differences, run:
```sh
git diff <file>
```

or if you'd wish to use a diff tool instead such as `meld`, type the following command just once

```sh
git config --global diff.tool meld
git config --global --add difftool.prompt false
```

Then, inspect a file using
```sh
git difftool <file>
```

or all changes in the root directory with `meld`
```sh
meld $CASSIOPEE
```

2. You can then add files and/or folders
```sh
git add <file1> <file2> <folder1> <file3>
```

or add all modified files/folders to your future commit:
```sh
git add -u
```

You are now ready to commit to your local repo - provide a meaningful description of the changes:
```sh
git commit -m "Module name: message"
```

3. Pushing will udpate your remote fork  
```sh
git push
```

<br> </br>

## Submitting changes to the upstream repository

Finally, to propose your modifications to the `dev` branch of the upstream/official repository, you have to 
create a Pull Request. Go to the webpage of your remote fork: https://github.com/github_username/Cassiopee/, click on the _Contribute_ button and follow the instructions (give a title to your PR, describe the modifications and review your changes).
