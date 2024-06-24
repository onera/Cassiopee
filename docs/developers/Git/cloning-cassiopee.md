# Cloning Cassiopee on Linux

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

5. Clone Cassiopee
```sh 
git clone git@github.com:onera/Cassiopee.git
```
