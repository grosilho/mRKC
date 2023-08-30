FROM alpine
USER root
WORKDIR /Libraries
RUN apk add --no-cache eigen-dev git build-base cmake make
RUN git clone https://gitlab.onelab.info/gmsh/gmsh.git
ADD https://www.fftw.org/fftw-3.3.10.tar.gz .
RUN chmod +x fftw-3.3.10.tar.gz && tar -xf fftw-3.3.10.tar.gz
RUN cd gmsh && mkdir build && cd build && \
    cmake -DENABLE_BUILD_DYNAMIC=1 .. && make && make install
RUN cd fftw-3.3.10 && mkdir build && cd build && cmake -S .. -B . && \
    make && make install 
WORKDIR /multirate
COPY . .
RUN mkdir build && cd build && cmake -S .. -B . && make
WORKDIR /multirate/build
USER root

