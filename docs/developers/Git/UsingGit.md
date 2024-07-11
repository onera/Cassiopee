# Using git

## Forks
To help developping Cassiopee, you have to first fork the
Cassiopee repository.

Regularly, you will have to "sync" your fork using the "sync" button
on the web site of your fork.
When syncing, your fork is updated from the main repository.

## Pulling
Pulling is updating your sources from your fork (after syncing):
```sh
git pull
```
If you have modified the same sources as the one that has been touched on the main repository, this operation will fail. You then have to move your modification aside by doing:
```sh
git stash
git pull
```
The pull will be ok, since no local modification exists. You now have to put back your modifications with:
```sh
git stash apply
```
You local modifications will be merged. Merge problem may occur and you will have to edit some files to choose the right code.

## Commiting
Commiting will commit your modifications locally. Nothing is sent. 
To check the files you have modified:
```sh
git status
```
To see the differences:
```sh
git diff <file>
```

Before commiting, you must first tag files for commit using:
```sh
git add <file>
```
or tag all modified files:
```sh
git add -u
```

You can then commit:
```sh
git commit -m "message"
```

## Pushing to your fork
Pushing is submitting your modifications to your fork. 
```sh
git push
```

## Merging to main repository

Finally, to send your modifications to the main repository, you have to 
submit a merge request. Go to your fork web site and click on "submit merge request". Add some comments about what you have down and click ok.

