FROM alpine
USER root
WORKDIR /Libraries
RUN apk add --no-cache eigen-dev git build-base cmake make
WORKDIR /multirate
COPY . .
RUN mkdir build && cd build && cmake -S .. -B . && make
WORKDIR /multirate/build
USER root
ENTRYPOINT ["/multirate/build/MultirateIntegrators"]
