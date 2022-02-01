FROM ubuntu:18.04

LABEL maintainer="sw@upmem.com"

ARG DEBIAN_FRONTEND=noninteractive
RUN adduser --disabled-password --home /home/docker-build --uid 10000 --gecos ",,," builder
RUN echo "builder ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers
RUN apt -y update
RUN apt -y install vim
RUN apt -y install apt-utils

# various utils
RUN apt -y install net-tools
RUN apt -y install dvipng texlive-latex-extra texlive ncftp
RUN apt -y install chrpath graphviz
RUN apt -y install doxygen
RUN apt -y install curl
RUN apt -y install wget
RUN apt install -y ninja-build libncurses-dev autoconf gperf bison flex texinfo help2man gawk libtool-bin git git-gui libconfig++-dev libz-dev valgrind cmake-curses-gui libboost-thread-dev libboost-system-dev 

# basic python on system
RUN apt -y install python lsb-release unzip
RUN apt -y install python-pip
RUN apt -y install python-dev
RUN apt -y install python-setuptools
RUN apt -y install virtualenv

# required python3.7 sytem alternative install 
RUN apt -y install libreadline-gplv2-dev libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libgdbm-dev libc6-dev libbz2-dev
RUN apt -y install libffi-dev
RUN apt -y install liblzma-dev lzma
RUN wget https://www.python.org/ftp/python/2.7/Python-2.7.tgz && tar -xzf Python-2.7.tgz && cd Python-2.7 && ./configure && make &&  make altinstall
RUN python2 -m pip install --user Cython

# nympy dev mode
RUN git clone https://github.com/numpy/numpy && cd numpy && git checkout v1.13.3 && git submodule update --init --recursive  && python setup.py install

# -------------------------------------------------------------------
# Create upmem user
RUN adduser --disabled-password upmem
RUN adduser upmem sudo
RUN echo 'upmem:upmem' | chpasswd
RUN echo 'upmem ALL=(ALL) NOPASSWD: ALL' >> /etc/sudoers

#  --------------------------------------------------------------------
#  Install upmem SDK not using docker cache to make sure wget is performed
ARG CACHEBUST=1
RUN mkdir -p /home/upmem
RUN wget -O /home/upmem/upmem_sdk_amd64.deb http://sdk-releases.upmem.com/2021.3.0/ubuntu_18.04/upmem_2021.3.0_amd64.deb 
RUN apt-get install -y /home/upmem/upmem_sdk_amd64.deb && \
    rm -f /home/upmem/upmem_sdk_amd64.deb

# cmake 3.21.4
RUN pip install cmake==3.21.4

# install specific perf module
RUN apt install -y python3-pip
RUN pip3 install pyyaml 'pyelftools>=0.26'
RUN git clone --depth 1 --branch "v$(uname -r | sed 's/\([0-9]\+\.[0-9]\+\).*/\1/g')" https://git.kernel.org/pub/scm/linux/kernel/git/torvalds/linux.git
RUN cd linux/tools/perf && make install && cp ./perf /usr/bin/perf
