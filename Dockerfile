FROM ubuntu:22.04

RUN apt-get update && \
    apt-get install -y build-essential wget bash cmake git iproute2 python3 python3-pip net-tools && \
    pip install tcconfig && \
    rm -rf /var/lib/apt/lists/*

#RUN apt-get remove -y cmake

WORKDIR /home

COPY ./install-dependencies-in-container.sh /home/install-dependencies-in-container.sh

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

#RUN apt-get update && \
#    apt-get install -y libboost-all-dev

#RUN wget https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.gz && \
#    tar -xzf boost_1_86_0.tar.gz && \
#    cd boost_1_86_0 && \
#    ./bootstrap.sh && \
#    ./b2 install && \
#    cd .. && \
#    rm -rf boost_1_86_0 boost_1_86_0.tar.gz