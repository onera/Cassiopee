# Running Cassiopee in a Docker container


## 1. Installing Docker and its dependencies

- Dependencies:
```sh
sudo apt install apt-transport-https ca-certificates curl software-properties-common
```

- Docker:
```sh
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt update && apt install docker-ce -y
sudo systemctl status docker
```

- Configure Docker:

```sh
sudo usermod -aG docker ${USER}
groups ${USER}
```

<br></br>

## 2. Pulling from DockerHub

- Please pull the official image of Cassiopee from DockerHub 

```sh
docker pull cassiopee486/cassiopee:<tag>
```

where _tag_ is a version tag such as `v4.0b` (in the remainder of these notes, `v4.0b` is used).
Available tags, **from `v4.0b` onwards**, are recalled [on this page](https://github.com/onera/Cassiopee/tags).

- Verify that the image is now present in the list of downloaded images

```sh
docker images
```

<br></br>

## 3. Running a container

- The Cassiopee container can be run interactively for version `v4.0b` like so

```sh
xhost +local:docker
docker run -it --rm --network=host --privileged --volume="$HOME/.Xauthority:/root/.Xauthority:rw" -v /tmp/.X11-unix:/tmp/.X11-unix -v /dev/dri:/dev/dri -e DISPLAY=unix$DISPLAY cassiopee486/cassiopee:v4.0b
```

and the instance of the container will be removed after it execution thanks to `--rm`.  
**Note that none of the modifications you may have made in the container will persist.** Please consider using _volumes_ or _bind mounts_ if this is something you may benefit from. An example is given below for a bind mount with read and write permissions using the command-line flag `-v /home/user/git/io:/io`, where the local folder `/home/user/git/io` is mapped to the `/io` folder in the container

```sh
xhost +local:docker
docker run -it --rm --network=host --privileged --volume="$HOME/.Xauthority:/root/.Xauthority:rw" -v /tmp/.X11-unix:/tmp/.X11-unix -v /dev/dri:/dev/dri -v /home/user/git/io:/io -e DISPLAY=unix$DISPLAY cassiopee486/cassiopee:v4.0b
```

- After execution, feel free to check the list of running instances (it should be none)
```sh
docker ps -a
```

<br></br>

## 4. Deleting a Docker image

To delete an outdated docker image, first list all existing images, copy the _hash_ of the image you would like to delete and remove from hash as

```sh
docker images
docker rmi <imageHash>
```
