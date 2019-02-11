FROM ubuntu:16.04
MAINTAINER olga.botvinnik@czbiohub.org

WORKDIR /tmp

USER root

# Install basics
ENV PACKAGES git make ca-certificates zlib1g-dev build-essential curl wget cmake apt-utils

### don't modify things below here for version updates etc.

WORKDIR /home

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

# Add add-apt-repository function
RUN apt-get update
RUN apt-get install -y software-properties-common

# Install gcc6 specifically
RUN add-apt-repository ppa:ubuntu-toolchain-r/test
RUN apt-get update && apt-get install -y g++-6
RUN g++ --version

# Install
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-6 60 --slave /usr/bin/g++ g++ /usr/bin/g++-6

WORKDIR /
RUN git clone https://github.com/dnbaker/dashing/
WORKDIR /dashing
RUN pwd
RUN make update dashing
RUN cp /dashing/dashing /bin

# Test that getting help on dashing command works
RUN dashing -h

WORKDIR /
