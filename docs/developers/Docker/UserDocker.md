# Running Cassiopee in a Docker container


1. Installing Docker and its dependencies

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

2. Pulling from DockerHub

- Please pull the official image of Cassiopee from DockerHub 

```sh
docker pull cassiopee/cassiopee:<tag>
```

where <tag> is a version tag such as `v4.0b`. In the remainder of these notes, `v4.0b` is used.
Available tags are recalled [on this page](https://github.com/onera/Cassiopee/tags), *from `v4.0b` onwards*.

- Verify that the image is now present in your list of downloaded images

```sh
docker images
```

3. Running a container

- The Cassiopee container can be run interactively like so

```sh
xhost +local:docker
docker run -it --rm -v /tmp/.X11-unix:/tmp/.X11-unix -v /dev/dri:/dev/dri -e DISPLAY=unix$DISPLAY cassiopee/cassiopee:main
```

and the instance of the container will be removed after it execution thanks to `--rm`.  
*Note that none of the modifications your may have made in the container will persist.* Please consider using _volumes_ or _bind mounts_ if this is something you may benefit from.

- After execution, feel free to check the list of running instances (it should be none).
```sh
docker ps -a
```

4. Deleting a Docker image

To delete an outdated docker image, first list all existing images, copy the _hash_ of the image you would like to delete and remove from hash as

```sh
docker images
docker rmi <imageHash>
```
