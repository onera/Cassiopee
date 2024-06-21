# -- Useful git shortcuts --

## Detailed & compact display of git log with graph
```shell
alias gitlog="git log --all --graph --pretty=format:'%C(auto)%d%Creset %C(Yellow) %h %Creset|| %C(Cyan) %as %Creset || %C(dim magenta) %an %Creset|| %s'"
```
See below for an example of the output of the above command:

![image](https://github.com/onera/Cassiopee/assets/124277807/fb1f2574-f793-4b11-afee-0cc563f10e17)

## Add branch detail to working path directory - in case of multiple branches
To be placed in the environment script & to be sourced
```shell
blackcolprompt='\[\033[30m\]'
purplecolprompt='\[\033[35m\]'
    
PS1='`if [ \$? = 0 ]; then echo \[\033[01\;32m\] ; else echo \[\033[01\;31m\]; fi`\u@\h\[\033[01;34m\] \w$(
  gitinfo="$(git branch 2>/dev/null | grep "^*" | cut -c3-)";
  if [ "$gitinfo" ]
  then
  	printf '$purplecolprompt'@"$gitinfo"'$redcolprompt'
  fi
  ) \$\[\033[00m\]'
```
The above command adds the '@main' seen below:

![image](https://github.com/onera/Cassiopee/assets/124277807/e15c130f-9b80-4a1c-b2b8-7c6a44e3584e)


## Git grep with line number
Git grep with line number for (1) a single variable or (2) a single variable in a specific working directory. The latter automatically ignores the */test directory and ignores the .rst files.

e.g. 
1) gitgrep var2
2) gitgrep var 2 Converter
```shell
#first argument: variable names
#second argument: local path to directory
gitgrep () {
    if [ "$#" -lt 2 ]
    then
	echo git grep -n "$1";
	git grep -n "$1";
    elif [ "$#" -lt 3 ]
    else
	echo git grep -n "$1" -- ':!'"$2"'/test/*'  ':!*.rst' "$2";
	git grep -n "$1" -- ':!'"$2"'/test/*'  ':!*.rst' "$2";
    else
	echo "gitgrep can only take a max of 2 arguments"
    fi
}
```
