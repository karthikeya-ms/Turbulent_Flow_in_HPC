# Docker

## Build Docker Image
```shell
docker build -t ns-eof .
```

## Run Docker image
```shell
docker run -it -v ${PWD}:/work --rm --privileged ns-eof /bin/bash
```

A prebuilt Docker image is available in [Dockerhub](https://hub.docker.com/r/tumi5/ns-eof).
