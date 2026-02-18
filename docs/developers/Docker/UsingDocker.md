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
docker pull cassiopee486/cassiopee:<tag>-<pyversion>
```

where _tag_ is a Cassiopee version tag (e.g.,`v4.1`) and _pyversion_ is a version of python (e.g., `py3.12.4`).
Available tags are recalled [on this page](https://hub.docker.com/r/cassiopee486/cassiopee/tags).
In the remainder of these notes, `v4.1-py3.12.4` is used.


- Verify that the image is now present in the list of downloaded images

```sh
docker images
```

<br></br>

## 3. Running a container

### 3.1 Manual configuration

The Cassiopee container can be run interactively for version `v4.1` with:

```sh
xhost +local:docker
docker run -it --rm --network=host --privileged \
    --volume="$HOME/.Xauthority:/root/.Xauthority:rw" \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    -v /dev/dri:/dev/dri \
    -e DISPLAY=unix$DISPLAY \
    cassiopee486/cassiopee:v4.1-py3.12.4
```

- The `--rm` flag ensures that the container is automatically removed after it exits.  
- Any changes made inside the container will not persist. To retain data, consider using volumes or bind mounts.  

For example, to map a local folder `/home/user/git/io` to `/io` in the container with read-write access:

```sh
xhost +local:docker
docker run -it --rm --network=host --privileged \
    --volume="$HOME/.Xauthority:/root/.Xauthority:rw" \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    -v /dev/dri:/dev/dri \
    -v $HOME/io:/io \
    -e DISPLAY=unix$DISPLAY \
    cassiopee486/cassiopee:v4.1-py3.12.4
```

After execution, you can verify that no containers are still running:

```sh
docker ps -a
```

### 3.2 Docker Compose

A `docker-compose.yml` file is provided in the root directory of Cassiopee, preconfigured with the same settings as in section 3.1.

To run the container interactively using Docker Compose:

```sh
xhost +local:docker
cd $CASSIOPEE
cp .env.template .env  # configure any environment variables if needed
docker compose run --rm cassiopee
```

Docker Compose automatically sets up the container with the correct volumes, display, and privileges.

As with manual execution, `--rm` ensures the container is removed when you exit.

Use the `.env` file to configure the tag and python version used without modifying the Compose file directly.

<br></br>

## 4. Deleting a Docker image

To delete an outdated docker image, first list all existing images, copy the _hash_ of the image you would like to delete and remove from hash as

```sh
docker images
docker rmi <imageHash>
```
