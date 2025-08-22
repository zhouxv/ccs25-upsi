FROM ubuntu:22.04

RUN apt-get update
RUN apt-get install -y net-tools iproute2 python3 python3-pip
RUN pip install tcconfig
RUN apt-get install -y build-essential wget bash cmake git

WORKDIR /home

COPY ./install-dependencies-in-container.sh /home/install-dependencies-in-container.sh

# note: must download boost_1_86_0.tar.bz2 before, hash value: 1bed88e40401b2cb7a1f76d4bab499e352fa4d0c5f31c0dbae64e24d34d7513b
COPY ./boost_1_86_0.tar.bz2 /home
RUN chmod +x /home/install-dependencies-in-container.sh

RUN /home/install-dependencies-in-container.sh

COPY ./sparseComp /home/sparseComp
COPY ./tests /home/tests
COPY ./build-bench.sh /home/build-bench.sh
COPY ./build-rls.sh /home/build-rls.sh
COPY ./build-tests.sh /home/build-tests.sh
COPY CMakeLists.txt /home/CMakeLists.txt

RUN chmod +x ./*.sh && \
    ./build-bench.sh && \
    cp ./build/fuzzylinf_bench ./
