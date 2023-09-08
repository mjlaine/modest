# modest libraries in a container
#
# Build the image with
#
#   docker build --rm -t modest .
#
# Start interactive shell:
#
#   docker run --rm -it -h modest modest
#
# Map current directory to "work" at container
#   docker run --rm -it -h modest -v $(pwd):/home/modest/work  modest
#
# quick test
#
#   docker run --rm -it modest sh -c 'cp -r /opt/modest/boxo . && cd boxo && make run'

ARG UBUNTU_RELEASE=22.04
FROM ubuntu:$UBUNTU_RELEASE
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

MAINTAINER Marko Laine <marko.laine@fmi.fi>

ARG DEBUG=yes

RUN useradd -ms /bin/bash modest
RUN usermod -aG sudo modest

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    vim unzip make git gfortran libopenblas-dev sudo \
    apt-transport-https ca-certificates \
    && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

# fetch code from github and build libraries
RUN cd /opt/ \
    && git clone --recurse-submodules https://github.com/mjlaine/modest.git \
    && cd modest \
    && make install -e DEBUG=$DEBUG

USER modest
WORKDIR /home/modest
