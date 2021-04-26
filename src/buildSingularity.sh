#!/usr/bin/env bash

# This script documents how to build the singularity container
# from the Dockerfile

if [[ $(/usr/bin/id -u) -ne 0 ]]; then
    echo "Must run script with sudo or as root"
    exit
fi

# exit on errors
trap 'exit' ERR

# build docker container
docker build -t gvanno .

# docker registry server to host docker image locally
# do nothing if running, otherwise try to start registry or create registry
[ $(docker inspect -f '{{.State.Running}}' registry) == "true" ] \
  || docker container start registry \
  || docker run -d -p 5000:5000 --restart=always --name registry registry:2

# push image to local registry
docker tag gvanno localhost:5000/gvanno
docker push localhost:5000/gvanno

# build singularity image from the docker cached local gvanno image
sudo SINGULARITY_NOHTTPS=1 singularity build gvanno.sif docker-daemon://gvanno:latest

# stop registry server
docker container stop registry
docker container rm registry
