# Cloning Cassiopee

## Add an SSH key to your Github account

1. Start the SSH agent    
```sh
eval `ssh-agent -s`
```

2. Create a new SSH key (no passphrase, add a meaningful suffix to the key) and add it to the SSH agent
```sh
ssh-keygen -t ed25519 -C <firstname.name@example.com>
ssh-add <full path to the private key (id_edXXXXX)>
```

3. Copy the contents of the public key (id_edXXXXX.pub) to your online Github account
```sh    
cat <full path to the public key>
```

Go to github.com/\<username\>, click on your profile in the top right-hand corner, then  
Settings > Access > SSH and GPG keys > New SSH key  
and paste it in the text box.

4. Set up your github account in a terminal window
```sh 
git config --global user.name <username>
git config --global user.email <firstname.name@example.com>
```

## Clone or Fork Cassiopee

1. If you are not planning on making any modifications to the source code, then you can clone the `main` branch of the Cassiopee repository.
Change directory to where you'd like the Cassiopee Github repository to be and type:

```sh 
git clone --branch main --single-branch git@github.com:onera/Cassiopee.git
```

2. To contribute to the Github repository instead, please consider creating a fork first: on the [github repo](https://github.com/onera/Cassiopee), click on the _Fork_ button followed by _Create new fork_. Then, you can clone the `dev` branch of your Fork as:

```sh 
git clone git@github.com:<github_username>/Cassiopee.git
```

The alias `origin` that is defined on your machine is pointing to your remote fork repository. Let's add another shortname called `upstream` to store the URL of the upstream/official repository:

```sh 
cd Cassiopee
git remote -v
git remote add upstream git@github.com:onera/Cassiopee.git
git remote -v
```

The next time you will want to synchronise your fork from the command line, please type

```sh 
git status # make sure you have no local changes
git pull upstream dev
git push
```

If you have local changes that you would wish to commit and push after synchronising your fork and local / working directories, please stash them as explained [here](https://github.com/onera/Cassiopee/blob/main/docs/developers/Git/UsingGit.md#stashing).

## Install Cassiopee

To install Cassiopee, please visit one of these Installation pages:  
- [Windows](https://github.com/onera/Cassiopee/blob/main/docs/developers/Install/msys2.md)
- [Ubuntu](https://github.com/onera/Cassiopee/blob/main/docs/developers/Install/ubuntu.md)
