FROM vibsinglecellnf/scanpy:0.6.2

ENV DEBIAN_FRONTEND=noninteractive
RUN BUILDPKGS="git" && \
    apt-get update && \
    apt-get install -y --no-install-recommends apt-utils debconf locales && dpkg-reconfigure locales && \
    apt-get install -y --no-install-recommends $BUILDPKGS

RUN git clone https://github.com/tabdelaal/SpaGE.git

ENV PYTHONPATH $PYTHONPATH:${PWD}/SpaGE/
